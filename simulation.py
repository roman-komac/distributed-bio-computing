from field import Field
from population import *
import numpy as np
from multiprocessing import Pipe
import os
from collections import defaultdict
import time
import signal
import atexit

class Simulation():
    def __init__(self, size_x=10, size_y=10, step=0.1, start_time=0.0, end_time=1000.0, multicore=False):
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
        self._cell_type_time = defaultdict(float)
        self._cell_type_num = defaultdict(int)
        self._field_time = 0.0
        self._field_num = 0
        self._switch_to_parallel = False

        if multicore:
            # Register step function
            self.step = self.__before_step_parallel
            # Current process is not child
            self._is_child = False
            # List of spawned child processes
            self._child_processes = list()
            # Register termination handlers
            signal.signal(signal.SIGINT, self.__signal_handler)
            signal.signal(signal.SIGTERM, self.__signal_handler)
            signal.signal(signal.SIGHUP, self.__signal_handler)
            atexit.register(self.__exit_handler)
        else:
            # Register step function
            self.step = self.__step

    def __exit_handler(self):
        if not self._is_child:
            # Kill child processes
            print("Number of child processes: " + str(len(self._child_processes)))
            for child_process_pid in self._child_processes:
                print("Killing process " + str(child_process_pid))
                os.kill(child_process_pid, signal.SIGTERM)
                pid, status = os.waitpid(child_process_pid, 0)
                if status == 0:
                    print("Process successfully terminated")
                else:
                    print("Process " + str(child_process_pid) + " orphaned")
            
            self._child_processes.clear()

            for pipe in self._field_pipes:
                pipe[0].close()
                pipe[1].close()
            for pipe in self._population_pipes:
                pipe[0].close()
                pipe[1].close()
            for pipe in self._pipes.values():
                pipe[0].close()
                pipe[1].close()
            # Clean up pipes and events

        exit()

    def __signal_handler(self, sig, frame):
        self.__exit_handler()

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
    
    # Huen method, single thread
    def __step(self):
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


    # Step per field, parallel
    def __step_field(self):
        while True:
            self._current_step = self._main_pipe.recv()

            curr_fields = {field_name : self._fields[field_name].get_current_state() for field_name in self._selected_keys}
            temp_fields = {field_name : self._fields[field_name].get_current_state() for field_name in self._selected_keys}

            for pipe in self._output_pipes:
                pipe.send([ (key, val["internal"]) for key,val in temp_fields.items() ])

            d_fields = {field_name : self._fields[field_name].step(self._current_time) for field_name in self._selected_keys}

            for pipe in self._input_pipes:
                pr = pipe.recv()
                for field_name, field_value in pr.items():
                    if field_name in d_fields:
                        d_fields[field_name]["internal"] += field_value

            for key in self._selected_keys:
                temp_fields[key]["internal"] += (self._step_size/2) * d_fields[key]["internal"]
                temp_fields[key]["external"] += (self._step_size/2) * d_fields[key]["external"]
                self._fields[key].set_current_state(temp_fields[key])

            for pipe in self._output_pipes:
                pipe.send([ (key, val["internal"]) for key,val in temp_fields.items() ])

            d_fields = {field_name : self._fields[field_name].step(self._current_time + self._step_size/2) for field_name in self._selected_keys}

            for pipe in self._input_pipes:
                for field_name, field_value in pipe.recv().items():
                    if field_name in d_fields:
                        d_fields[field_name]["internal"] += field_value
        
            for key in self._selected_keys:
                curr_fields[key]["internal"] += self._step_size * d_fields[key]["internal"]
                curr_fields[key]["external"] += self._step_size * d_fields[key]["external"]
                self._fields[key].set_current_state(curr_fields[key])

            self._current_step += 1
            self._current_time = self._current_step * self._step_size

            self._main_pipe.send(curr_fields)


    # Step per population, parallel
    def __step_population(self):
        while True:
            self._current_step = self._main_pipe.recv()

            curr_populations = {population_name : self._populations[population_name].get_current_state() for population_name in self._selected_keys}
            temp_populations = {population_name : self._populations[population_name].get_current_state() for population_name in self._selected_keys}

            for pipe in self._input_pipes:
                pr = pipe.recv()
                for field_name, field_value in pr:
                    self._fields[field_name]._internal = field_value 

            d_populations = {population_name : self._populations[population_name].step(self._current_time) for population_name in self._selected_keys}
            
            d_fields = dict()
            for d_population in d_populations.values():
                if "output_fields" in d_population.keys():
                    out_fields = list(d_population["output_fields"].items())
                    for (out_field_name, out_field) in out_fields:
                        if out_field_name not in d_fields.keys():
                            d_fields[out_field_name] = out_field
                            continue
                        d_fields[out_field_name] += out_field
            
            for pipe in self._output_pipes:
                pipe.send(d_fields)

            for key in self._selected_keys:
                for state_key in temp_populations[key]["mRNA"].keys():
                    temp_populations[key]["mRNA"][state_key] += (self._step_size/2) * d_populations[key]["mRNA"][state_key]
                for state_key in temp_populations[key]["protein"].keys():
                    temp_populations[key]["protein"][state_key] += (self._step_size/2) * d_populations[key]["protein"][state_key]
                self._populations[key].set_current_state(temp_populations[key])

            
            for pipe in self._input_pipes:
                for field_name, field_value in pipe.recv():
                    self._fields[field_name]._internal = field_value 

            d_populations = {population_name : self._populations[population_name].step(self._current_time + self._step_size/2) for population_name in self._selected_keys}

            d_fields = dict()
            for d_population in d_populations.values():
                if "output_fields" in d_population.keys():
                    out_fields = list(d_population["output_fields"].items())
                    for (out_field_name, out_field) in out_fields:
                        if out_field_name not in d_fields.keys():
                            d_fields[out_field_name] = out_field
                            continue
                        d_fields[out_field_name] += out_field

            for pipe in self._output_pipes:
                pipe.send(d_fields)

            for key in self._selected_keys:
                for state_key in curr_populations[key]["mRNA"].keys():
                    curr_populations[key]["mRNA"][state_key] += self._step_size * d_populations[key]["mRNA"][state_key]
                for state_key in curr_populations[key]["protein"].keys():
                    curr_populations[key]["protein"][state_key] += self._step_size * d_populations[key]["protein"][state_key]
                self._populations[key].set_current_state(curr_populations[key])

            
            self._main_pipe.send(curr_populations)


    # Prepare processes to switch to parallel implementation
    def __prepare_for_parallel(self):
        cell_type_timing = dict()
        for key in self._cell_type_time.keys():
            cell_type_timing[key] = self._cell_type_time[key] / self._cell_type_num[key]

        lst = [(pop_name,cell_type_timing[pop._cell_type]) for pop_name,pop in self._populations.items()]
        sorted_lst = sorted(lst, key=lambda x: x[1], reverse=True)

        sm_fields = len(self._fields) * self._field_time / self._field_num
        a = defaultdict(float)
        for k,v in self._populations.items():
            a[v._cell_type] += cell_type_timing[v._cell_type]
        sm_populations = sum(a.values())

        # Fields take longer to process
        if sm_fields > sm_populations and round(sm_fields/sm_populations) > 1:
            coefficient = min(round(sm_fields/sm_populations), 3)
            self._field_pipes = [Pipe() for v in range(coefficient)]
            self._field_keys = []
            cf = len(self._fields)/coefficient
            m = list()
            fkeys = list(self._fields.keys())
            for i in range(coefficient):
                s = round(i*cf)
                e = round((i+1)*cf)
                self._field_keys.append(set(fkeys[s:e]))

            self._population_pipes = [Pipe()]
            self._population_keys = [set(self._populations.keys())]

        # Populations take longer to process
        elif round(sm_populations/sm_fields) > 1:
            coefficient = min(round(sm_populations/sm_fields), 3)
            self._field_pipes = [Pipe()]
            self._field_keys = [set(self._fields.keys())]
            self._population_pipes = [Pipe() for v in range(coefficient)]
            self._population_keys = [set() for c in range(coefficient)]
            # balance
            rsk = [0.0 for i in range(coefficient)]
            l = 0
            for pop_name,pop_time in sorted_lst:
                rsk[l%coefficient] += pop_time
                self._population_keys[l%coefficient].add(pop_name)
                if rsk[(l-1)%coefficient] < rsk[l%coefficient]:
                    l += 1

        # Both take approx the same time
        else:
            self._field_pipes = [Pipe()]
            self._field_keys = [set(self._fields.keys())]
            self._population_pipes = [Pipe()]
            self._population_keys = [set(self._populations.keys())]

        # Initialize pipes for communication between field/population processes
        self._pipes = dict()
        for field_i in range(len(self._field_pipes)):
            for population_i in range(len(self._population_pipes)):
                self._pipes[(field_i,population_i)] = Pipe()

        # For every field division fork child process
        for field_i in range(len(self._field_pipes)):
            field_i_pipes = [self._pipes[(field_i,r)] for r in range(len(self._population_pipes))]
            self._input_pipes  = [field_i_pipe[0] for field_i_pipe in field_i_pipes]
            self._output_pipes = [field_i_pipe[0] for field_i_pipe in field_i_pipes]
            self._main_pipe = self._field_pipes[field_i][0]
            self._selected_keys = self._field_keys[field_i]

            pid = os.fork()
            if pid == 0:
                self._is_child = True
                self.__step_field()
            else:
                self._child_processes.append(pid)
                continue

        # For every population division fork child process
        for population_i in range(len(self._population_pipes)):
            population_i_pipes = [self._pipes[(r,population_i)] for r in range(len(self._field_pipes))]
            self._input_pipes  = [population_i_pipe[1] for population_i_pipe in population_i_pipes]
            self._output_pipes = [population_i_pipe[1] for population_i_pipe in population_i_pipes]
            self._main_pipe = self._population_pipes[population_i][0]
            self._selected_keys = self._population_keys[population_i]

            pid = os.fork()
            if pid == 0:
                self._is_child = True
                self.__step_population()
            else:
                self._child_processes.append(pid)
                continue


    # Parallel method step
    def __step_parallel(self):
        for pipe in self._field_pipes:
            pipe[1].send(self._current_step)

        for pipe in self._population_pipes:
            pipe[1].send(self._current_step)
        
        for pipe in self._field_pipes:
            for key,value in pipe[1].recv().items():
                self._fields[key].set_current_state(value)
        for pipe in self._population_pipes:
            for key,value in pipe[1].recv().items():
                self._populations[key].set_current_state(value)

        self._current_step += 1
        self._current_time = self._current_step * self._step_size

        return self._current_time >= self._end_time


    # This method is used to estimate time per population/field
    def __before_step_parallel(self):
        if self._switch_to_parallel:
            return self.__step_parallel()

        current_populations_timing = { key : self._cell_type_time[key] / max(self._cell_type_num[key], 1) for key in self._cell_type_num.keys() }
        current_field_timing = self._field_time / max(self._field_num, 1)

        curr_populations = {population_name : population.get_current_state() for (population_name, population) in self._populations.items()}
        curr_fields = {field_name : field.get_current_state() for (field_name, field) in self._fields.items()}

        d_populations = dict() 
        for (population_name, population) in self._populations.items():
            pt = time.time()
            d_populations[population_name] = population.step(self._current_time)
            self._cell_type_num[population._cell_type] += 1
            self._cell_type_time[population._cell_type] += time.time() - pt

        d_fields = dict()
        for (field_name, field) in self._fields.items():
            pt = time.time()
            d_fields[field_name] = field.step(self._current_time)
            self._field_num += 1
            self._field_time += time.time() - pt
        
        
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


        d_populations = dict() 
        for (population_name, population) in self._populations.items():
            pt = time.process_time()
            d_populations[population_name] = population.step(self._current_time + self._step_size/2)
            self._cell_type_num[population._cell_type] += 1
            self._cell_type_time[population._cell_type] += time.process_time() - pt

        d_fields = dict()
        for (field_name, field) in self._fields.items():
            pt = time.process_time()
            d_fields[field_name] = field.step(self._current_time + self._step_size/2)
            self._field_num += 1
            self._field_time += time.process_time() - pt

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

        diff_populations_timing = { key : self._cell_type_time[key] / max(self._cell_type_num[key], 1) for key in self._cell_type_num.keys() }
        diff_field_timing = self._field_time / max(self._field_num, 1)

        mn = - 0.1 ** 9
        mx = 0.1 ** 9

        if all( mn < (cpt - diff_populations_timing[key]) < mx for key, cpt in current_populations_timing.items() ) and (mn < (current_field_timing - diff_field_timing) < mx):
            self._switch_to_parallel = True
            self.__prepare_for_parallel()

        return self._current_time >= self._end_time
