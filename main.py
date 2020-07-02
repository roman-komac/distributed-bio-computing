#!/usr/bin/env python3

# math library
import math

# matplotlib plotter
import matplotlib.pyplot as plt

# pyqtgraph plotter
import pyqtgraph as pg

# PyLab
from pylab import *

# Simulation
from simulation import Simulation

# Populations
from population import *

# Fields
from field import Field

# Timing library
import time



# Create simulation with fields of size 10x10, step 0.1 min
# start_time defaults to 0 min, end_time equals 200 hours
sim = Simulation(size_x=10, size_y=10, step=0.1, end_time=60 * 500)

# Field dict
fields = dict()

# Clock field
fields["clk"] = Field("clk")

# FLIP-FLOP 1 Q FIELD
fields["ff1"] = Field("ff1")
# FLIP-FLOP 2 Qc FIELD
fields["ff2"] = Field("ff2")

# Decoder 1 output FIELD
fields["dec1"] = Field("dec1")
# Decoder 2 output FIELD
fields["dec2"] = Field("dec2")
# Decoder 3 output FIELD
fields["dec3"] = Field("dec3")
# Decoder 4 output FIELD
fields["dec4"] = Field("dec4")


# Population dict
populations = dict()

# Create REPRESSILATOR
populations["rep"] = Repressilator("rep")
populations["rep"].num_cells = 10
populations["rep"].alpha *= 4
populations["rep"].beta *= 0.75
populations["rep"].kappa = 0.05
populations["rep"].add_input_field("CLK", fields["clk"])
populations["rep"].add_output_field("A", fields["clk"])
populations["rep"].init_random_scale = 30
populations["rep"].init_offset_scale = 40

# Create FLIP-FLOP 1
populations["ff1"] = FlipFlop("ff1")
populations["ff1"].num_cells = 15
populations["ff1"].add_input_field("CLK", fields["clk"])
populations["ff1"].add_input_field("D", fields["ff2"])
populations["ff1"].add_output_field("Q", fields["ff1"])
populations["ff1"].init_random_scale = 20
populations["ff1"].init_offset_scale = 50

# Create FLIP-FLOP 2
populations["ff2"] = FlipFlop("ff2")
populations["ff2"].num_cells = 15
populations["ff2"].add_input_field("CLK", fields["clk"])
populations["ff2"].add_input_field("D", fields["ff1"])
populations["ff2"].add_output_field("Qc", fields["ff2"])
populations["ff2"].init_random_scale = 20
populations["ff2"].init_offset_scale = 50


# Create OR 1 (Decoder 1)
populations["dec1"] = Or("dec1")
populations["dec1"].num_cells = 5
populations["dec1"].add_input_field("D1", fields["ff1"])
populations["dec1"].add_input_field("D2", fields["ff2"])
populations["dec1"].add_input_field("SYNC", fields["dec1"])
populations["dec1"].add_output_field("A", fields["dec1"])

# Create single inverse OR (Decoder 2)
populations["dec2"] = SInvOr("dec2")
populations["dec2"].num_cells = 5
populations["dec2"].add_input_field("D1", fields["ff1"])
populations["dec2"].add_input_field("D2", fields["ff2"])
populations["dec2"].add_input_field("SYNC", fields["dec2"])
populations["dec2"].add_output_field("A", fields["dec2"])

# Create single inverse OR (Decoder 3)
populations["dec3"] = SInvOr("dec3")
populations["dec3"].num_cells = 5
populations["dec3"].add_input_field("D1", fields["ff2"])
populations["dec3"].add_input_field("D2", fields["ff1"])
populations["dec3"].add_input_field("SYNC", fields["dec3"])
populations["dec3"].add_output_field("A", fields["dec3"])

# Create double inverse OR (Decoder 4)
populations["dec4"] = DInvOr("dec4")
populations["dec4"].num_cells = 5
populations["dec4"].add_input_field("D1", fields["ff1"])
populations["dec4"].add_input_field("D2", fields["ff2"])
populations["dec4"].add_input_field("SYNC", fields["dec4"])
populations["dec4"].add_output_field("A", fields["dec4"])


# Create signal inverser (Decoder 1)
populations["dec_not1"] = Not("dec_not1")
populations["dec_not1"].num_cells = 1
populations["dec_not1"].add_input_field("D", fields["dec1"])

# Create signal inverser (Decoder 1)
populations["dec_not2"] = Not("dec_not2")
populations["dec_not2"].num_cells = 1
populations["dec_not2"].add_input_field("D", fields["dec2"])

# Create signal inverser (Decoder 1)
populations["dec_not3"] = Not("dec_not3")
populations["dec_not3"].num_cells = 1
populations["dec_not3"].add_input_field("D", fields["dec3"])

# Create signal inverser (Decoder 1)
populations["dec_not4"] = Not("dec_not4")
populations["dec_not4"].num_cells = 1
populations["dec_not4"].add_input_field("D", fields["dec4"])


# Add populations and fields to simulation
sim.add_populations(populations)
sim.add_fields(fields)


steps_per_record = 20
n_points = math.floor( (sim._end_time - sim._start_time) / sim._step_size) // steps_per_record

# Record population concentrations
record_full = dict()

for (population_name, protein) in [("rep","A"), ("ff1","Q"), ("ff1","Qc"), ("ff2","Q"), ("ff2","Qc")]:
    record_full[(population_name,protein)] = np.zeros((n_points, populations[population_name].num_cells))

for population_name in ["dec1", "dec2", "dec3", "dec4"]:
    record_full[(population_name,"A")] = np.zeros((n_points, populations[population_name].num_cells))

for population_name in ["dec_not1", "dec_not2", "dec_not3", "dec_not4"]:
    record_full[(population_name,"A")] = np.zeros((n_points, populations[population_name].num_cells))

int_dec1_full  = np.zeros((n_points, populations["dec_not1"].num_cells))
ext_dec1_full  = np.zeros((n_points, populations["dec_not1"].num_cells))

int_ff1_full  = np.zeros((n_points, populations["dec1"].num_cells))
ext_ff1_full  = np.zeros((n_points, populations["dec1"].num_cells))

tstart = time.time()

# Run step-by-step
# Returns false when out of steps
curr_hour = 0
while not sim.step():
    if sim._current_step % steps_per_record == 0:
        step10 = sim._current_step // steps_per_record
        
        for (population_name, protein) in record_full.keys():
            population_state = populations[population_name].get_current_state()["protein"][protein]
            record_full[(population_name, protein)][step10,:] = [ population_state[id] for id in sim._cell_ids[population_name] ]

        dec1_field_state = fields["dec1"].get_current_state()
        ff1_field_state  = fields["ff1"].get_current_state()
        int_dec1_full[step10,:] = [ dec1_field_state["internal"][id] for id in sim._cell_ids["dec_not1"] ]
        ext_dec1_full[step10,:] = [ dec1_field_state["external"][id] for id in sim._cell_ids["dec_not1"] ]
        int_ff1_full[step10,:]  = [ ff1_field_state["internal"][id] for id in sim._cell_ids["dec1"] ]
        ext_ff1_full[step10,:]  = [ ff1_field_state["external"][id] for id in sim._cell_ids["dec1"] ]

    if sim._current_time >= curr_hour * 60:
        curr_hour += 1
        sys.stdout.write('\033[2K\033[1G')
        print("  " + str(sim._current_time/60) + " h", end="\r")

print("Simulation took " + ("%0.1f" % (time.time() - tstart)) + " seconds")

fig,(row1,row2,row3,row4) = subplots(4,4)

T = np.linspace(0, (sim._end_time - sim._start_time)/60, n_points)
row1[0].plot(T, record_full[("rep","A")])
row1[0].set_title('concentrations of A_osc (repressilator)')

row2[0].plot(T, int_dec1_full)
row2[0].plot(T, ext_dec1_full)
row2[0].set_title('int concentrations of dec1 autoinducer')

row3[0].plot(T, int_ff1_full)
row3[0].plot(T, ext_ff1_full)
row3[0].set_title('int concentrations of ff1 autoinducer')

row1[1].plot(T, record_full[("ff1","Q")])
row1[1].set_title('concentrations of FF1 Q')

row2[1].plot(T, record_full[("ff1","Qc")])
row2[1].set_title('concentrations of FF1 Qc')

row3[1].plot(T, record_full[("ff2","Q")])
row3[1].set_title('concentrations of FF2 Q')

row4[1].plot(T, record_full[("ff2","Qc")])
row4[1].set_title('concentrations of FF2 Qc')


row1[2].plot(T, record_full[("dec1","A")])
row1[2].set_title('concentrations of Decoder 1 A')

row2[2].plot(T, record_full[("dec2","A")])
row2[2].set_title('concentrations of Decoder 2 A')

row3[2].plot(T, record_full[("dec3","A")])
row3[2].set_title('concentrations of Decoder 3 A')

row4[2].plot(T, record_full[("dec4","A")])
row4[2].set_title('concentrations of Decoder 4 A')


row1[3].plot(T, record_full[("dec_not1","A")])
row1[3].set_title('concentrations of NOT Decoder 1 A')

row2[3].plot(T, record_full[("dec_not2","A")])
row2[3].set_title('concentrations of NOT Decoder 2 A')

row3[3].plot(T, record_full[("dec_not3","A")])
row3[3].set_title('concentrations of NOT Decoder 3 A')

row4[3].plot(T, record_full[("dec_not4","A")])
row4[3].set_title('concentrations of NOT Decoder 4 A')

plt.show()

exit()
