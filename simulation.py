from field import Field
from population import *
import numpy as np
import concurrent.futures
from multiprocessing import Pool, Manager
import os
from collections import defaultdict

class Simulation():
    def __init__(self, size_x=10, size_y=10, step=0.1, start_time=0.0, end_time=1000.0):
        self._fields = dict()
        self._populations = dict()
        self._cell_positions = dict()
        self._cell_positions_all = None
        self._cell_ids = dict()
        self._num_cells_cumulative = 0
        self._step_size = step
        self._current_step = 0
        self._size = (size_y, size_x)
        self._cell_positions_all = np.zeros(self._size)
        self._start_time = self._current_time = start_time
        self._end_time = end_time
            

    def add_population(self, population_name, population):
        if population_name is None:
            raise RuntimeError('Population name not specified')
        
        if not isinstance(population, Population):
            raise RuntimeError('Parameter population is not an instance of Population')
        
        env_size = self._size[0] * self._size[1]
        if population.num_cells + self._num_cells_cumulative > env_size:
            raise RuntimeError('Number of cells of population exceeds the number of available spaces')
        
        self.__initialize_population(population_name, population)

    def remove_population(self, population_name):
        if population_name is None:
            raise RuntimeError('Population name not specified')

        if population_name not in self._populations.keys():
            raise RuntimeError('Population does not exist in environment')

        cell_pos = self._cell_positions.pop(population_name)
        self._cell_positions_all -= cell_pos
        self._cell_positions_all = self._cell_positions_all.clip(min=0)
        self._populations.pop(population_name)
        self._cell_ids.pop(population_name)
            
    def add_populations(self, populations):
        for (population_name, population) in populations.items():
            self.add_population(population_name, population)

    def add_field(self, field_name, field):
        if field_name is None:
            raise RuntimeError('Field name not specified')
        
        if not isinstance(field, Field):
            raise RuntimeError('Parameter field is not an instance of Field')
        
        self.__initialize_field(field_name, field)
        

    def remove_field(self, field_name):
        if field_name is None:
            raise RuntimeError('Field name not specified')

        if field_name not in self._fields.keys():
            raise RuntimeError('Field does not exist in environment')

        self._fields.pop(field_name)
            
    def add_fields(self, fields):
        for (field_name, field) in fields.items():
            self.add_field(field_name, field)

    def __initialize_population(self, population_name, population):
        self._populations[population_name] = population

        free_positions = np.nonzero(self._cell_positions_all == 0)
        free_positions_len = free_positions[0].size
        selected_positions = np.random.choice(free_positions_len, population.num_cells, replace=False)
        
        self._cell_ids[population_name] = [pos for pos in zip(free_positions[0][selected_positions],free_positions[1][selected_positions])]
        self._cell_positions[population_name] = np.zeros(self._size)
        for pos in self._cell_ids[population_name]:
            self._cell_positions_all[pos] = 1
            self._cell_positions[population_name][pos] = 1

        self._cell_ids[population_name].sort( key=lambda id: id[0]*self._size[0] + id[1] )
        self._populations[population_name].initialize_states(self._size, self._cell_positions[population_name])
        
    def __initialize_field(self, field_name, field):
        self._fields[field_name] = field
        self._fields[field_name].initialize_states(self._size, self._cell_positions_all)
    
    def step(self):

        # Huen method
        curr_populations = {population_name : population.get_current_state() for (population_name, population) in self._populations.items()}
        curr_fields = {field_name : field.get_current_state() for (field_name, field) in self._fields.items()}
        
        d_populations = {population_name : population.step(self._current_time) for (population_name, population) in self._populations.items()}
        d_fields = {field_name : field.step(self._current_time) for (field_name, field) in self._fields.items()}
        for d_population in d_populations.values():
            if "output_fields" in d_population.keys():
                out_fields = list(d_population["output_fields"].items())
                for (out_field_name, out_field) in out_fields:
                    d_fields[out_field_name]["internal"] += out_field
        

        temp_populations = {population_name : population.get_current_state() for (population_name, population) in self._populations.items()}
        temp_fields = {field_name : field.get_current_state() for (field_name, field) in self._fields.items()}

        for key in self._populations.keys():
            for state_key in temp_populations[key]["mRNA"].keys():
                temp_populations[key]["mRNA"][state_key] += (self._step_size/2) * d_populations[key]["mRNA"][state_key]
            for state_key in temp_populations[key]["protein"].keys():
                temp_populations[key]["protein"][state_key] += (self._step_size/2) * d_populations[key]["protein"][state_key]
            self._populations[key].set_current_state(temp_populations[key])

        for key in self._fields.keys():
            temp_fields[key]["internal"] += (self._step_size/2) * d_fields[key]["internal"]
            temp_fields[key]["external"] += (self._step_size/2) * d_fields[key]["external"]
            self._fields[key].set_current_state(temp_fields[key])


        d_populations = {population_name : population.step(self._current_time + self._step_size/2) for (population_name, population) in self._populations.items()}
        d_fields = {field_name : field.step(self._current_time + self._step_size/2) for (field_name, field) in self._fields.items()}
        for d_population in d_populations.values():
            if "output_fields" in d_population.keys():
                out_fields = list(d_population["output_fields"].items())
                for (out_field_name, out_field) in out_fields:
                    d_fields[out_field_name]["internal"] += out_field


        for key in self._populations.keys():
            for state_key in curr_populations[key]["mRNA"].keys():
                curr_populations[key]["mRNA"][state_key] += self._step_size * d_populations[key]["mRNA"][state_key]
            for state_key in curr_populations[key]["protein"].keys():
                curr_populations[key]["protein"][state_key] += self._step_size * d_populations[key]["protein"][state_key]
            self._populations[key].set_current_state(curr_populations[key])

        for key in self._fields.keys():
            curr_fields[key]["internal"] += self._step_size * d_fields[key]["internal"]
            curr_fields[key]["external"] += self._step_size * d_fields[key]["external"]
            self._fields[key].set_current_state(curr_fields[key])

        self._current_step += 1
        self._current_time = self._current_step * self._step_size

        return self._current_time >= self._end_time
