from file_handling import HDF5File, XDMFCheckpointHandler
from optimizer import ScipyOptimizer
from fields import DDDbFields
from space_definition import Spaces
from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from fields.field_updates import FieldUpdates
import numpy as np



class DDDbTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 file: fenics.XDMFFile,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: DDDbFields,
                 mesh: fenics.Mesh,
                 spaces: Spaces,
                 strain_file_name: str,
                 material_parameters_file_name: str,
                 initial_material_parameters: np.ndarray,
                 checkpoint_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields)
        self.checkpoint_file = XDMFCheckpointHandler(file_name=checkpoint_file_name, append_to_existing=False,
                                                     field=fields.imported_displacement_field,
                                                     field_name=fields.u_new.name())

        self.tensor_space = spaces.tensor_space
        # self.t2d = fenics.vertex_to_dof_map(spaces.tensor_space)
        self.number_of_vertices = mesh.num_vertices()
        self.strain_file_name = strain_file_name
        self.material_parameters_file_name = material_parameters_file_name
        self.fields = fields
        initial_values_for_optimizer = np.random.random(
            self.fields.new_constitutive_relation_multiplicative_parameters.function_space().dim())
        self.optimizer = ScipyOptimizer(fields=self.fields, initial_values=initial_values_for_optimizer,
                                        fem_solver=fem_solver, spaces=spaces)



    def transform_material_parameters(self):
        strain_vec = fenics.project(fenics.sym(fenics.grad(self.fields.u_new)), V=self.tensor_space).vector()[:]
        strain_tens = np.reshape(strain_vec, newshape=(self.number_of_vertices, 2, 2))
        parameters_vec = fenics.project(self.fields.new_constitutive_relation_multiplicative_parameters,
                                        V=self.tensor_space).vector()[:]
        parameters_tens = np.reshape(parameters_vec, newshape=(self.number_of_vertices, 2, 2))
        return strain_tens, parameters_tens


    def run(self, i: int):
        print('iteration: {}'.format(i))
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t_float, i)
        self.checkpoint_file.load(i)
        self.optimizer.run()
        self.field_updates.run(fields=self.fields)
        self.file.write(self.fields.u_new, (i + 1)*self.time_params.delta_t_float)



    def close(self):
        self.checkpoint_file.close()