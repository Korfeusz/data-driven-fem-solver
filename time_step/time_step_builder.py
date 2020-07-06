from .time_step import TimeStep
from .dddb_time_step import DDDbTimeStep
from .elastodynamics_time_step import ElastodynamicsTimeStep
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
        self.mesh = None
        self.hdf5_file_name = None

    def set(self,  alpha_params=None, time_params=None, fem_solver=None,
                          file=None, boundary_excitation=None, field_updates=None, fields=None,
            mesh=None, hdf5_file_name=None):
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
        if mesh is not None:
            self.mesh = mesh
        if hdf5_file_name is not None:
            self.hdf5_file_name = hdf5_file_name


    def build(self):
        if self.time_step_type == ElastodynamicsTimeStep:
            return ElastodynamicsTimeStep(self.alpha_params, self.time_params, self.fem_solver,
                                          self.file, self.boundary_excitation, self.field_updates,
                                          self.fields, self.mesh, self.hdf5_file_name)
        if self.time_step_type == DDDbTimeStep:
            return DDDbTimeStep(self.alpha_params, self.time_params, self.fem_solver,
                                self.file, self.boundary_excitation, self.field_updates, self.fields,
                                self.mesh, self.hdf5_file_name)