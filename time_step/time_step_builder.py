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
        self.checkpoint_file_name = None
        self.spaces = None
        self.strain_file_name = None
        self.material_parameters_file_name = None
        self.initial_material_parameters = None

    def set(self, alpha_params=None, time_params=None, fem_solver=None,
            file=None, boundary_excitation=None, field_updates=None, fields=None,
            mesh=None, checkpoint_file_name=None, spaces=None, strain_file_name=None, material_parameters_file_name=None,
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
        if checkpoint_file_name is not None:
            self.checkpoint_file_name = checkpoint_file_name
        if spaces is not None:
            self.spaces = spaces
        if strain_file_name is not None:
            self.strain_file_name = strain_file_name
        if material_parameters_file_name is not None:
            self.material_parameters_file_name = material_parameters_file_name
        if initial_material_parameters is not None:
            self.initial_material_parameters = initial_material_parameters


    def build(self):
        if self.time_step_type == ElastodynamicsTimeStep:
            return ElastodynamicsTimeStep(self.alpha_params, self.time_params, self.fem_solver,
                                          self.file, self.boundary_excitation, self.field_updates,
                                          self.fields, checkpoint_file_name=self.checkpoint_file_name)
        if self.time_step_type == DDDbTimeStep:
            return DDDbTimeStep(alpha_params=self.alpha_params, time_params=self.time_params, fem_solver=self.fem_solver,
                                file=self.file, boundary_excitation=self.boundary_excitation,
                                field_updates=self.field_updates, fields=self.fields,
                                checkpoint_file_name=self.checkpoint_file_name, mesh=self.mesh, spaces=self.spaces,
                                strain_file_name=self.strain_file_name,
                                material_parameters_file_name=self.material_parameters_file_name,
                                initial_material_parameters=self.initial_material_parameters)