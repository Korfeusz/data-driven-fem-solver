import fenics
import numpy as np


class TimeSteppingParameters:
    def __init__(self, total_time, number_of_steps):
        self.total_time = total_time
        self.number_of_steps = number_of_steps

    @property
    def delta_t(self):
        return fenics.Constant(self.total_time/self.number_of_steps)

    @property
    def linear_time_space(self):
        return np.linspace(0, self.total_time, self.number_of_steps+1)

    @property
    def delta_t_float(self):
        return float(self.delta_t)