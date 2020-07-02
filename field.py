# Numpy
import numpy as np

class Field():
    def __init__(self, name):
        if name is not None:
            self.name = name
        else:
            raise RuntimeError("Field name cannot be None")

        # Parameters
        self.field_diff_rate = 0.5
        self.membrane_diff_rate = 1
        self.kSi = 1 # Degradation in mol/min
        self.kSe = 0.01 # Degradation in mol/min
        self._h2 = 0.25
        self._size = None
        self.border_fixed = False
        self.periodic_bounds = True
        self._internal = None
        self._external = None

    def initialize_states(self, size, cell_positions):
        self._size = size
        self._cell_positions = cell_positions
        self._external = np.random.random_sample(self._size)
        self._internal = np.zeros(self._size)
    
    def get_current_state(self):
        return {
            "internal" : self._internal, 
            "external" : self._external
        }

    def set_current_state(self, states):
        self._internal = states["internal"]
        self._external = states["external"]

    def step(self, t):
        membrane_diff = self.membrane_diff_rate * (self._internal - self._external)
        d_internal = - self.kSi * self._internal - self._cell_positions * membrane_diff
        d_external = - self.kSe * self._external + self._cell_positions * membrane_diff

        if self.periodic_bounds:
            # Periodic bounds wrap-around
            external_xx = self.field_diff_rate * ( np.roll(self._external, 1, 1) + np.roll(self._external, -1, 1)  -2 * self._external )/self._h2
            external_yy = self.field_diff_rate * ( np.roll(self._external, 1, 0) + np.roll(self._external, -1, 0)  -2 * self._external )/self._h2
        else:
            # Create padded matrix to incorporate Neumann boundary conditions 
            size_y = self._size[0]
            size_x = self._size[1]
            p_external = np.zeros((size_y + 2, size_x + 2))
            p_external[1:(size_y + 1), 1:(size_x + 1)] = self._external
            p_external[0, :] = p_external[2, :]
            p_external[(size_y + 1), :] = p_external[(size_y - 1), :]
            p_external[:, 0] = p_external[:, 2]
            p_external[:, (size_x + 1)] = p_external[:, (size_x - 1)]

            # Calculate diffusion part of the equations
            external_xx = self.field_diff_rate * (p_external[1:(size + 1), 0:size] + p_external[1:(size + 1), 2:(size + 2)] -2 * self._external)/self._h2
            external_yy = self.field_diff_rate * (p_external[0:size, 1:(size + 1)] + p_external[2:(size + 2), 1:(size + 1)] -2 * self._external)/self._h2

        d_external += external_xx + external_yy

        if self.border_fixed:
            # Border as distortion centers
            d_external[1:size,(size-1)] = 0
            d_external[1:size,0] = 0
            d_external[(size-1),1:size] = 0
            d_external[0,1:size] = 0

        return {"internal":d_internal , "external":d_external}