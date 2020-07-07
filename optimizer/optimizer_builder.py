from typing import Type

from .optimizer import Optimizer
from .scipy_optimizer import ScipyOptimizer


class OptimizerBuilder:
    def __init__(self, optimizer_type: Type[Optimizer]):
        self.optimizer_type = optimizer_type
        self.alpha_params = None
        self.time_params = None
        self.fem_solver = None
        self.file = None
        self.boundary_excitation = None
        self.field_updates = None
        self.fields = None
        self.mesh = None
        self.hdf5_file_name = None
        self.spaces = None
        self.strain_file_name = None
        self.material_parameters_file_name = None
        self.initial_material_parameters = None

    def set(self,  alpha_params=None, time_params=None, fem_solver=None,
                          file=None, boundary_excitation=None, field_updates=None, fields=None,
            mesh=None, hdf5_file_name=None, spaces=None, strain_file_name=None, material_parameters_file_name=None,
            initial_material_parameters=None):
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
        if self.spaces is not None:
            self.spaces = spaces
        if self.strain_file_name is not None:
            self.strain_file_name = strain_file_name
        if self.material_parameters_file_name is not None:
            self.material_parameters_file_name = material_parameters_file_name
        if self.initial_material_parameters is not None:
            self.initial_material_parameters = initial_material_parameters


    def build(self):
        if self.optimizer_type == ScipyOptimizer:
            return ScipyOptimizer()