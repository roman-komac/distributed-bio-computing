# Numpy
import numpy as np

class Population:
    def __init__(self, cell_type=None, name=None):
        if cell_type is not None:
            self._cell_type = cell_type
        else:
            raise RuntimeError("Cell type cannot be None")

        if name is not None:
            self.name = name
        else:
            raise RuntimeError("Population name cannot be None")
    
        self._mRNA    = dict()
        self._protein = dict()

        self._input_fields  = dict()
        self._output_fields = dict()

        self._cell_positions = None
        self.num_cells = 0
        
        self._size = (10,10)
        self.init_random_scale = 100
        self.init_offset_scale = 10

    def add_input_field(self, protein_name, field):
        if protein_name is None:
            raise RuntimeError('Protein name not specified')
        
        if field is None:
            raise RuntimeError('Input field not specified')

        self._input_fields[protein_name] = field

    def add_output_field(self, protein_name, field):
        if protein_name is None:
            raise RuntimeError('Output field name not specified')
        
        if field is None:
            raise RuntimeError('Output field not specified')

        self._output_fields[protein_name] = field

    def remove_input_field(self, field):
        if field is None:
            raise RuntimeError('Input field not specified')
        
        self._input_fields.pop(field.name)

    def remove_output_field(self, field):
        if field is None:
            raise RuntimeError('Output field not specified')
        
        self._output_fields.pop(field.name)

    def clear_input_fields(self):
        self._input_fields = dict()

    def clear_output_fields(self):
        self._output_fields = dict()

    def get_current_state(self):
        return {
            "mRNA"    : self._mRNA, 
            "protein" : self._protein
        }

    def set_current_state(self, states):
        self._mRNA = states["mRNA"]
        self._protein = states["protein"]


class Repressilator(Population):
    def __init__(self, name):
        super().__init__(cell_type='Repressilator', name=name)
        # Parameters
        self.Kd      = 10.0   # threshold
        self.n       = 2      # exponent
        self.alpha   = 1.25   # mol/min
        self.alpha0  = 0.0025 # mol/min
        self.beta    = 0.1    # mol/min
        self.delta_p = 0.01   # mol/min
        self.delta_m = 0.01   # mol/min
        self.kappa   = 0.02   # mol/min
        self.kS      = 0.5    # mol/min

    def initialize_states(self, size, cell_positions):
        self._cell_positions = cell_positions
        self._size = size
        self._mRNA    = { mRNA : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale) for mRNA in ["mA", "mB", "mC"]}
        self._protein = { protein : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale) for protein in ["A", "B", "C"]}

    def step(self, t):
        CLK_state = self._input_fields["CLK"].get_current_state()
        CLK = self._cell_positions * CLK_state["internal"]

        # mA
        d_mA = self._cell_positions * (self.alpha/(1 + (self._protein["C"]/self.Kd) ** self.n) + self.alpha0 - self.delta_m * self._mRNA["mA"])
        # mB
        d_mB = self._cell_positions * (self.alpha/(1 + (self._protein["A"]/self.Kd) ** self.n) + self.alpha0 - self.delta_m * self._mRNA["mB"])
        # mC
        d_mC = self._cell_positions * (self.alpha/(1 + (self._protein["B"]/self.Kd) ** self.n) + self.alpha0 - self.delta_m * self._mRNA["mC"] + (self.kappa * CLK)/(1 + CLK))
        # A
        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"] - self.delta_p * self._protein["A"])
        # B
        d_B  = self._cell_positions * (self.beta * self._mRNA["mB"] - self.delta_p * self._protein["B"])
        # C
        d_C  = self._cell_positions * (self.beta * self._mRNA["mC"] - self.delta_p * self._protein["C"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA, "mB":d_mB, "mC":d_mC}, 
            "protein" : {"A" :d_A,  "B" :d_B,  "C" :d_C},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }


class FlipFlop(Population):
    def __init__(self, name):
        super().__init__(cell_type='FlipFlop', name=name)
        # Parameters
        self.Kd      = 1.04   # threshold
        self.n       = 4.35   # exponent
        self.alpha_1 = 0.0781 # mol/min
        self.alpha_2 = 0.1851 # mol/min
        self.alpha_3 = 0.0736 # mol/min
        self.alpha_4 = 0.1858 # mol/min
        self.delta_1 = 0.0322 # mol/min
        self.delta_2 = 0.0115 # mol/min
        self.beta    = 0.02   # mol/min
        self.delta_m = 0.01   # mol/min
        self.kS      = 1.1    # mol/min

    def initialize_states(self, size, cell_positions):
        self._cell_positions = cell_positions
        self._size = size
        self._mRNA    = { mRNA : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale) for mRNA in ["mA", "mAc", "mQ", "mQc"]}
        self._protein = { protein : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale) for protein in ["A", "Ac", "Q", "Qc"]}

    def step(self, t):
        CLK_state = self._input_fields["CLK"].get_current_state()
        CLK = self._cell_positions * CLK_state["internal"]

        D_state = self._input_fields["D"].get_current_state()
        D = self._cell_positions * D_state["internal"]

        #A
        d_mA  =  self.alpha_1 * ( ((D/self.Kd) ** self.n)/(1 + ((D/self.Kd) ** self.n) + ((CLK/self.Kd) ** self.n) + ((D/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n)))
        d_mA  += self.alpha_2 * (1/(1 + ((self._protein["Ac"]/self.Kd) ** self.n))) - self.delta_1 * self._mRNA["mA"]
        #Ac
        d_mAc =  self.alpha_1 * (1/(1 + ((D/self.Kd) ** self.n) + ((CLK/self.Kd) ** self.n) + ((D/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n))) 
        d_mAc += self.alpha_2 * (1/(1 + ((self._protein["A"]/self.Kd) ** self.n)))  - self.delta_1 * self._mRNA["mAc"]
        #Q
        d_mQ  =  self.alpha_3 * ((((self._protein["A"]/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n))/(1 + ((self._protein["A"]/self.Kd) ** self.n) + ((CLK/self.Kd) ** self.n) + ((self._protein["A"]/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n))) 
        d_mQ  += self.alpha_4 * (1/(1 + ((self._protein["Qc"]/self.Kd) ** self.n))) - self.delta_2 * self._mRNA["mQ"]
        #Qc
        d_mQc =  self.alpha_3 * ((((self._protein["Ac"]/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n))/(1 + ((self._protein["Ac"]/self.Kd) ** self.n) + ((CLK/self.Kd) ** self.n) + ((self._protein["Ac"]/self.Kd) ** self.n)*((CLK/self.Kd) ** self.n))) 
        d_mQc += self.alpha_4 * (1/(1 + ((self._protein["Q"]/self.Kd) ** self.n))) - self.delta_2 * self._mRNA["mQc"]

        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"]  - self.delta_m * self._protein["A"])
        d_Ac = self._cell_positions * (self.beta * self._mRNA["mAc"] - self.delta_m * self._protein["Ac"])
        d_Q  = self._cell_positions * (self.beta * self._mRNA["mQ"]  - self.delta_m * self._protein["Q"])
        d_Qc = self._cell_positions * (self.beta * self._mRNA["mQc"] - self.delta_m * self._protein["Qc"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA, "mAc":d_mAc, "mQ":d_mQ, "mQc":d_mQc},
            "protein" : {"A" :d_A,  "Ac" :d_Ac,  "Q" :d_Q,  "Qc" :d_Qc},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }


class Or(Population):
    def __init__(self, name, cell_type='Or'):
        super().__init__(cell_type=cell_type, name=name)
        # Parameters
        self.Kd      = 1.5   # threshold
        self.n       = 4      # exponent
        self.alpha   = 0.3    # mol/min
        self.beta    = 0.0125 # mol/min
        self.delta_p = 0.01   # mol/min
        self.delta_m = 0.01   # mol/min
        self.kS      = 1.5    # mol/min
        self.kappa   = 0.1    # mol/min

    def initialize_states(self, size, cell_positions):
        self._cell_positions = cell_positions
        self._size = size
        self._mRNA    = { "mA" : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale)}
        self._protein = { "A"  : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale)}

    def step(self, t):
        D1_state = self._input_fields["D1"].get_current_state()
        D1 = self._cell_positions * D1_state["internal"]

        D2_state = self._input_fields["D2"].get_current_state()
        D2 = self._cell_positions * D2_state["internal"]

        SYNC_state = self._input_fields["SYNC"].get_current_state()
        SYNC = self._cell_positions * SYNC_state["internal"]

        #mA
        d_mA = self._cell_positions * self.alpha * ( ((D1/self.Kd) ** self.n)/(1 + ((D1/self.Kd) ** self.n)) + ((D2/self.Kd) ** self.n)/(1 + ((D2/self.Kd) ** self.n)) ) + ((self.kappa * SYNC)/(1 + SYNC)) - self.delta_m * self._mRNA["mA"]
        # A
        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"] - self.delta_p * self._protein["A"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA},
            "protein" : {"A" :d_A},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }


class SInvOr(Or):
    def __init__(self, name):
        super().__init__(cell_type='SInvOr', name=name)
        self.alpha   = 0.15   # mol/min
        self.beta    = 0.025  # mol/min

    def step(self, t):
        D1_state = self._input_fields["D1"].get_current_state()
        D1 = self._cell_positions * D1_state["internal"]

        D2_state = self._input_fields["D2"].get_current_state()
        D2 = self._cell_positions * D2_state["internal"]

        SYNC_state = self._input_fields["SYNC"].get_current_state()
        SYNC = self._cell_positions * SYNC_state["internal"]

        #mA
        d_mA = self._cell_positions * self.alpha * ( 1/(1 + ((D1/self.Kd) ** self.n)) + ((D2/self.Kd) ** self.n)/(1 + ((D2/self.Kd) ** self.n)) ) + ((self.kappa)/(1 + SYNC)) - self.delta_m * self._mRNA["mA"]
        # A
        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"] - self.delta_p * self._protein["A"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA},
            "protein" : {"A" :d_A},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }

class DInvOr(Or):
    def __init__(self, name):
        super().__init__(cell_type='DInvOr', name=name)

    def step(self, t):
        D1_state = self._input_fields["D1"].get_current_state()
        D1 = self._cell_positions * D1_state["internal"]

        D2_state = self._input_fields["D2"].get_current_state()
        D2 = self._cell_positions * D2_state["internal"]

        SYNC_state = self._input_fields["SYNC"].get_current_state()
        SYNC = self._cell_positions * SYNC_state["internal"]

        #mA
        d_mA = self._cell_positions * self.alpha * ( 1/(1 + ((D1/self.Kd) ** self.n)) + 1/(1 + ((D2/self.Kd) ** self.n)) ) + ((self.kappa * SYNC)/(1 + SYNC)) - self.delta_m * self._mRNA["mA"]
        # A
        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"] - self.delta_p * self._protein["A"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA},
            "protein" : {"A" :d_A},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }


class Not(Population):
    def __init__(self, name):
        super().__init__(cell_type="Not", name=name)
        # Parameters
        self.Kd      = 0.75   # threshold
        self.n       = 4      # exponent
        self.alpha   = 1.5    # mol/min
        self.beta    = 0.25   # mol/min
        self.delta_p = 0.05   # mol/min
        self.delta_m = 0.05   # mol/min
        self.kS      = 0.5    # mol/min
        self.kappa   = 0.25   # mol/min

    def initialize_states(self, size, cell_positions):
        self._cell_positions = cell_positions
        self._size = size
        self._mRNA    = { "mA" : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale)}
        self._protein = { "A"  : self._cell_positions * (np.random.random_sample(self._size) * self.init_random_scale + self.init_offset_scale)}

    def step(self, t):
        D_state = self._input_fields["D"].get_current_state()
        D = self._cell_positions * D_state["internal"]

        SYNC_state = self._input_fields["SYNC"].get_current_state()
        SYNC = self._cell_positions * SYNC_state["internal"]

        #mA
        d_mA = self._cell_positions * ( self.alpha/(1 + ((D/self.Kd) ** self.n)) ) + ((self.kappa * SYNC)/(1 + SYNC)) - self.delta_m * self._mRNA["mA"]
        # A
        d_A  = self._cell_positions * (self.beta * self._mRNA["mA"] - self.delta_p * self._protein["A"])

        # Get output field
        out_fields = list(self._output_fields.items())
        return {
            "mRNA"    : {"mA":d_mA},
            "protein" : {"A" :d_A},
            "output_fields" : {f.name : (self._protein[p] * self.kS) for (p,f) in out_fields }
        }