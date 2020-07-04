import abc
from typing import Type


class TimeStep(abc.ABC):
    @abc.abstractmethod
    def __init__(self, alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields):
        self.alpha_params = alpha_params
        self.time_params = time_params
        self.fem_solver = fem_solver
        self.file = file
        self.boundary_excitation = boundary_excitation
        self.field_updates = field_updates
        self.fields = fields

    @abc.abstractmethod
    def run(self, i):
        pass
