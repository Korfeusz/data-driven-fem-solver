from file_handling import HDF5File, XDMFCheckpointHandler, NPYFile
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
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: DDDbFields,
                 mesh: fenics.Mesh,
                 spaces: Spaces,
                 dddb_output_file: str,
                 initial_material_parameters: np.ndarray,
                 in_checkpoint_file_name: str,
                 out_checkpoint_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)
        self.in_checkpoint_file = XDMFCheckpointHandler(file_name=in_checkpoint_file_name, append_to_existing=False,
                                                        field=fields.imported_displacement_field,
                                                        field_name=fields.u_new.name())
        self.out_checkpoint_file = XDMFCheckpointHandler(file_name=out_checkpoint_file_name, append_to_existing=False,
                                                        field=fields.u_new,
                                                        field_name=fields.u_new.name())
        self.v_out_checkpoint_file = XDMFCheckpointHandler(file_name='v_{}'.format(out_checkpoint_file_name), append_to_existing=False,
                                                        field=fields.v_old,
                                                        field_name=fields.v_old.name())
        self.a_out_checkpoint_file = XDMFCheckpointHandler(file_name='a_{}'.format(out_checkpoint_file_name), append_to_existing=False,
                                                        field=fields.a_old,
                                                        field_name=fields.a_old.name())
        self.np_files = NPYFile(file_name=dddb_output_file)
        self.tensor_space = spaces.tensor_space
        self.number_of_vertices = mesh.num_vertices()
        self.fields = fields
        initial_values_for_optimizer = np.random.random(
            self.fields.new_constitutive_relation_multiplicative_parameters.function_space().dim())
        self.optimizer = ScipyOptimizer(fields=self.fields, initial_values=initial_values_for_optimizer,
                                        fem_solver=fem_solver, spaces=spaces)



    def save_strain_and_params(self, iteration: int) -> None:
        strain_vec = fenics.project(fenics.sym(fenics.grad(self.fields.u_new)), V=self.tensor_space).vector()[:]
        strain_tens = np.reshape(strain_vec, newshape=(self.number_of_vertices, 2, 2))
        parameters = fenics.project(self.fields.new_constitutive_relation_multiplicative_parameters,
                                        V=self.fields.new_constitutive_relation_multiplicative_parameters.function_space()).vector()[:]
        if self.fields.new_constitutive_relation_multiplicative_parameters.function_space().dim() == self.fields.tensor_space.dim():
            parameters = np.reshape(parameters, newshape=(self.number_of_vertices, 2, 2))
        self.np_files.write(prefix='strain', iteration=iteration, array=strain_tens)
        self.np_files.write(prefix='params', iteration=iteration, array=parameters)


    def run(self, i: int):
        print('iteration: {}'.format(i))
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t_float, i)
        self.in_checkpoint_file.load(i)
        self.optimizer.save_file = XDMFCheckpointHandler(file_name='optimizer_save_file_{}.xdmf'.format(i), append_to_existing=False,
                                                         field=self.fields.u_new, field_name=self.fields.u_new.name())
        self.optimizer.run()
        self.v_out_checkpoint_file.write(i)
        self.a_out_checkpoint_file.write(i)
        self.save_strain_and_params(i)
        self.field_updates.run(fields=self.fields)
        self.out_checkpoint_file.write(i)
        if not self.optimizer.keep_going:
            self.halt = True



    def close(self):
        self.in_checkpoint_file.close()
        self.out_checkpoint_file.close()
        self.a_out_checkpoint_file.close()
        self.v_out_checkpoint_file.close()