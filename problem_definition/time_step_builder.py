from .time_step import TimeStep
from . import ElastodynamicsTimeStep
from typing import Type

class TimeStepBuilder:
    def __init__(self, time_step_type: Type[TimeStep]):
        self.time_step_type = time_step_type
        self.alpha_params = None
        self.time_params = None
        self.fem_solver = None
        self.file = None
        self.boundary_excitation = None
        self.field_updates = None
        self.fields = None

    def set(self,  alpha_params=None, time_params=None, fem_solver=None,
                          file=None, boundary_excitation=None, field_updates=None, fields=None):
        if alpha_params is not None:
            self.alpha_params = alpha_params
        if time_params is not None:
            self.time_params = time_params
        if fem_solver is not None:
            self.fem_solver = fem_solver
        if file is not None:
            self.file = file
        if boundary_excitation is not None:
            self.boundary_excitation = boundary_excitation
        if field_updates is not None:
            self.field_updates = field_updates
        if fields is not None:
            self.fields = fields

    def build(self):
        if self.time_step_type == ElastodynamicsTimeStep:
            return ElastodynamicsTimeStep(self.alpha_params, self.time_params, self.fem_solver,
                                          self.file, self.boundary_excitation, self.field_updates, self.fields)