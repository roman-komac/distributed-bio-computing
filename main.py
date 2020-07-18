#!/usr/bin/env python3

# math library
from math import floor

# matplotlib plotter
import matplotlib.pyplot as plt

# pyqtgraph plotter
import pyqtgraph as pg

# import system specific
import sys

# PyLab
from pylab import subplots

# Timing library
import time

# Simulation
from simulation import Simulation

# Populations
from population import *

# Fields
from field import Field


# Create simulation with fields of size 10x10, step 0.1 min
# start_time defaults to 0 min, end_time equals 1500 hours
# if available, the parallel version will be used to speed-up simulation (multicore)
sim = Simulation(size_x=10, size_y=10, step=0.1, end_time=60 * 1500, multicore=True)

# Dictionary of fields
fields = dict()

# Clock field (population rep A protein controlled autoinducer)
fields["clk"] = Field("clk")

# Flip-flop 1 field (population ff1 Q protein controlled autoinducer)
fields["ff1"] = Field("ff1")
# Flip-flop 2 field (population ff1 Qc protein controlled autoinducer)
fields["ff2"] = Field("ff2")

# Decoder 1 output field (population dec1 A protein controlled autoinducer)
fields["dec1"] = Field("dec1")
# Decoder 2 output field (population dec2 A protein controlled autoinducer)
fields["dec2"] = Field("dec2")
# Decoder 3 output field (population dec3 A protein controlled autoinducer)
fields["dec3"] = Field("dec3")
# Decoder 4 output field (population dec4 A protein controlled autoinducer)
fields["dec4"] = Field("dec4")


# Negator 1 output field (population dec_not1 A protein controlled autoinducer)
fields["dec_not1"] = Field("dec_not1")
# Negator 2 output field (population dec_not2 A protein controlled autoinducer)
fields["dec_not2"] = Field("dec_not2")
# Negator 3 output field (population dec_not3 A protein controlled autoinducer)
fields["dec_not3"] = Field("dec_not3")
# Negator 4 output field (population dec_not4 A protein controlled autoinducer)
fields["dec_not4"] = Field("dec_not4")


# Dictionary of populations
populations = dict()

# Create Repressilator population
populations["rep"] = Repressilator("rep")
populations["rep"].num_cells = 20
populations["rep"].kS *= 0.333
populations["rep"].add_input_field("CLK", fields["clk"])
populations["rep"].add_output_field("A", fields["clk"])
populations["rep"].init_random_scale = 70
populations["rep"].init_offset_scale = 20

# Create Flip-Flop 1 population
populations["ff1"] = FlipFlop("ff1")
populations["ff1"].num_cells = 15
populations["ff1"].add_input_field("CLK", fields["clk"])
populations["ff1"].add_input_field("D", fields["ff2"])
populations["ff1"].add_output_field("Q", fields["ff1"])
populations["ff1"].init_random_scale = 40
populations["ff1"].init_offset_scale = 10

# Create Flip-Flop 2 population
populations["ff2"] = FlipFlop("ff2")
populations["ff2"].num_cells = 15
populations["ff2"].add_input_field("CLK", fields["clk"])
populations["ff2"].add_input_field("D", fields["ff1"])
populations["ff2"].add_output_field("Qc", fields["ff2"])
populations["ff2"].init_random_scale = 40
populations["ff2"].init_offset_scale = 10


# Create Or population ( ff1 OR ff2 )
populations["dec1"] = Or("dec1")
populations["dec1"].num_cells = 5
populations["dec1"].Kd *= 1.25
populations["dec1"].add_input_field("D1", fields["ff1"])
populations["dec1"].add_input_field("D2", fields["ff2"])
populations["dec1"].add_input_field("SYNC", fields["dec1"])
populations["dec1"].add_output_field("A", fields["dec1"])
populations["dec1"].init_random_scale = 10
populations["dec1"].init_offset_scale = 10

# Create single inverted Or population ( (NOT ff1) OR ff2 )
populations["dec2"] = SInvOr("dec2")
populations["dec2"].num_cells = 5
populations["dec2"].add_input_field("D1", fields["ff1"])
populations["dec2"].add_input_field("D2", fields["ff2"])
populations["dec2"].add_input_field("SYNC", fields["dec2"])
populations["dec2"].add_output_field("A", fields["dec2"])
populations["dec2"].init_random_scale = 40
populations["dec2"].init_offset_scale = 10

# Create single inverted Or population ( ff1 OR (NOT ff2) )
populations["dec3"] = SInvOr("dec3")
populations["dec3"].num_cells = 5
populations["dec3"].add_input_field("D1", fields["ff2"])
populations["dec3"].add_input_field("D2", fields["ff1"])
populations["dec3"].add_input_field("SYNC", fields["dec3"])
populations["dec3"].add_output_field("A", fields["dec3"])
populations["dec3"].init_random_scale = 40
populations["dec3"].init_offset_scale = 10

# Create double inverted Or population ( (NOT ff1) OR (NOT ff2) )
populations["dec4"] = DInvOr("dec4")
populations["dec4"].num_cells = 5
populations["dec4"].Kd *= 0.6
populations["dec4"].add_input_field("D1", fields["ff1"])
populations["dec4"].add_input_field("D2", fields["ff2"])
populations["dec4"].add_input_field("SYNC", fields["dec4"])
populations["dec4"].add_output_field("A", fields["dec4"])
populations["dec4"].init_random_scale = 40
populations["dec4"].init_offset_scale = 10


# Create signal inverter population ( NOT dec1 )
populations["dec_not1"] = Not("dec_not1")
populations["dec_not1"].num_cells = 7
populations["dec_not1"].Kd = 0.7
populations["dec_not1"].add_input_field("D", fields["dec1"])
populations["dec_not1"].add_input_field("SYNC", fields["dec_not1"])
populations["dec_not1"].add_output_field("A", fields["dec_not1"])
populations["dec_not1"].init_random_scale = 40
populations["dec_not1"].init_offset_scale = 10

# Create signal inverter population ( NOT dec2 )
populations["dec_not2"] = Not("dec_not2")
populations["dec_not2"].num_cells = 7
populations["dec_not2"].Kd = 0.55
populations["dec_not2"].add_input_field("D", fields["dec2"])
populations["dec_not2"].add_input_field("SYNC", fields["dec_not2"])
populations["dec_not2"].add_output_field("A", fields["dec_not2"])
populations["dec_not2"].init_random_scale = 40
populations["dec_not2"].init_offset_scale = 10

# Create signal inverter population ( NOT dec3 )
populations["dec_not3"] = Not("dec_not3")
populations["dec_not3"].num_cells = 7
populations["dec_not3"].Kd = 0.55
populations["dec_not3"].add_input_field("D", fields["dec3"])
populations["dec_not3"].add_input_field("SYNC", fields["dec_not3"])
populations["dec_not3"].add_output_field("A", fields["dec_not3"])
populations["dec_not3"].init_random_scale = 40
populations["dec_not3"].init_offset_scale = 10

# Create signal inverter population ( NOT dec4 )
populations["dec_not4"] = Not("dec_not4")
populations["dec_not4"].num_cells = 7
populations["dec_not4"].Kd = 0.85
populations["dec_not4"].add_input_field("D", fields["dec4"])
populations["dec_not4"].add_input_field("SYNC", fields["dec_not4"])
populations["dec_not4"].add_output_field("A", fields["dec_not4"])
populations["dec_not4"].init_random_scale = 40
populations["dec_not4"].init_offset_scale = 10

# Add populations and fields to simulation
sim.add_populations(populations)
sim.add_fields(fields)


# Record population concentrations
record_populations_full = dict()

# Record field concentrations on the positions of population cells
record_fields_full = dict()

# Number of steps before the concentration is recorded
steps_per_record = 50

# Number of points to record
n_points = int(floor( (sim._end_time - sim._start_time) / sim._step_size) // steps_per_record)

# Populations and which proteins to record
for (population_name, protein) in [("rep","A"), ("ff1","Q"), ("ff1","Qc"), ("ff2","Q"), ("ff2","Qc")]:
    record_populations_full[(population_name,protein)] = np.zeros((n_points, populations[population_name].num_cells))

for population_name in ["dec1", "dec2", "dec3", "dec4", "dec_not1", "dec_not2", "dec_not3", "dec_not4"]:
    record_populations_full[(population_name,"A")] = np.zeros((n_points, populations[population_name].num_cells))

# Populations and which fields to record
for (population_name, field_name) in [("dec1", "ff1"),("dec1", "ff2"),("ff1", "clk"),("ff2", "clk")]:
    record_fields_full[(population_name, field_name)] = {
        "internal" : np.zeros((n_points, populations[population_name].num_cells)), 
        "external" : np.zeros((n_points, populations[population_name].num_cells))
    }


# Start timing simulation
tstart = time.time()

# Run step-by-step
# Call to sim.step() returns false when out of steps
curr_hour = 0
while not sim.step():
    if sim._current_step % steps_per_record == 0:
        # Step to record
        record_step = sim._current_step // steps_per_record
        
        # Record concentrations of proteins for populations
        for (population_name, protein) in record_populations_full.keys():
            population_state = populations[population_name].get_current_state()["protein"][protein]
            record_populations_full[(population_name, protein)][record_step,:] = [ population_state[id] for id in sim._cell_ids[population_name] ]

        # Record concentrations of fields at the positions of population cells
        for (population_name, field_name) in record_fields_full.keys():
            field_state = fields[field_name].get_current_state()
            record_fields_full[(population_name, field_name)]["internal"][record_step,:] = [ field_state["internal"][id] for id in sim._cell_ids[population_name] ]
            record_fields_full[(population_name, field_name)]["external"][record_step,:] = [ field_state["external"][id] for id in sim._cell_ids[population_name] ]

    # Every 30 minutes of simulation print current time
    if sim._current_time >= curr_hour:
        curr_hour += 30 #minutes
        sys.stdout.write('\033[2K\033[1G') # clear line and go to beginning
        print("%0.1f" % (sim._current_time/60) + " h", end="\r")


print("Simulation took " + ("%0.1f" % (time.time() - tstart)) + " seconds")

# Plot results
fig,(row1,row2,row3,row4) = subplots(4,4)

T = np.linspace(0, (sim._end_time - sim._start_time)/60, n_points)
row1[0].plot(T, record_populations_full[("rep","A")])
row1[0].set_title('concentrations of A protein in rep population')

row2[0].plot(T, record_fields_full[("ff1","clk")]["internal"])
row2[0].plot(T, record_fields_full[("ff1","clk")]["external"])
row2[0].set_title('int and ext concentrations of clk autoinducer in ff1 population')

row3[0].plot(T, record_fields_full[("ff2","clk")]["internal"])
row3[0].plot(T, record_fields_full[("ff2","clk")]["external"])
row3[0].set_title('int and ext concentrations of clk autoinducer in ff2 population')


row1[1].plot(T, record_populations_full[("ff1","Q")])
row1[1].set_title('concentrations of ff1 Q')

row2[1].plot(T, record_populations_full[("ff1","Qc")])
row2[1].set_title('concentrations of ff1 Qc')

row3[1].plot(T, record_populations_full[("ff2","Q")])
row3[1].set_title('concentrations of ff2 Q')

row4[1].plot(T, record_populations_full[("ff2","Qc")])
row4[1].set_title('concentrations of ff2 Qc')


row1[2].plot(T, record_populations_full[("dec1","A")])
row1[2].set_title('concentrations of dec1 A')

row2[2].plot(T, record_populations_full[("dec2","A")])
row2[2].set_title('concentrations of dec2 A')

row3[2].plot(T, record_populations_full[("dec3","A")])
row3[2].set_title('concentrations of dec3 A')

row4[2].plot(T, record_populations_full[("dec4","A")])
row4[2].set_title('concentrations of dec4 A')


row1[3].plot(T, record_populations_full[("dec_not1","A")])
row1[3].set_title('concentrations of dec_not1 A')

row2[3].plot(T, record_populations_full[("dec_not2","A")])
row2[3].set_title('concentrations of dec_not2 A')

row3[3].plot(T, record_populations_full[("dec_not3","A")])
row3[3].set_title('concentrations of dec_not3 A')

row4[3].plot(T, record_populations_full[("dec_not4","A")])
row4[3].set_title('concentrations of dec_not4 A')

plt.show()
