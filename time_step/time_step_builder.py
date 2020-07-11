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
        self.boundary_excitation = None
        self.field_updates = None
        self.fields = None
        self.mesh = None
        self.checkpoint_file_name = None
        self.in_checkpoint_file_name = None
        self.out_checkpoint_file_name = None
        self.spaces = None
        self.dddb_output_file_name = None
        self.initial_material_parameters = None

    def set(self, alpha_params=None, time_params=None, fem_solver=None,
            boundary_excitation=None, field_updates=None, fields=None,
            mesh=None, checkpoint_file_name=None, spaces=None, dddb_output_file=None,
            initial_material_parameters=None, in_checkpoint_file_name=None, out_checkpoint_file_name=None):
        if alpha_params is not None:
            self.alpha_params = alpha_params
        if time_params is not None:
            self.time_params = time_params
        if fem_solver is not None:
            self.fem_solver = fem_solver
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
        if dddb_output_file is not None:
            self.dddb_output_file_name = dddb_output_file
        if initial_material_parameters is not None:
            self.initial_material_parameters = initial_material_parameters
        if in_checkpoint_file_name is not None:
            self.in_checkpoint_file_name = in_checkpoint_file_name
        if out_checkpoint_file_name is not None:
            self.out_checkpoint_file_name = out_checkpoint_file_name



    def build(self):
        if self.time_step_type == ElastodynamicsTimeStep:
            return ElastodynamicsTimeStep(alpha_params=self.alpha_params, time_params=self.time_params,
                                          fem_solver=self.fem_solver,
                                          boundary_excitation=self.boundary_excitation, field_updates=self.field_updates,
                                          fields=self.fields, checkpoint_file_name=self.checkpoint_file_name)
        if self.time_step_type == DDDbTimeStep:
            return DDDbTimeStep(alpha_params=self.alpha_params, time_params=self.time_params, fem_solver=self.fem_solver,
                                boundary_excitation=self.boundary_excitation,
                                field_updates=self.field_updates, fields=self.fields,
                                in_checkpoint_file_name=self.in_checkpoint_file_name,
                                out_checkpoint_file_name=self.out_checkpoint_file_name,
                                mesh=self.mesh, spaces=self.spaces,
                                dddb_output_file=self.dddb_output_file_name,
                                initial_material_parameters=self.initial_material_parameters
                                )