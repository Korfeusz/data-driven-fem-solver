import fenics

from fem_solver import FemSolver
from fields import FieldUpdates, DDDbFields
from file_handling import XDMFCheckpointHandler
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation
from time_step import TimeStep
from time_stepping_parameters import TimeSteppingParameters
from data_driven_solver import DDDbLoader, ParameterUpdates
import numpy as np


class DDSolverTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: DDDbFields,
                 out_checkpoint_file_name: str,
                 dddb_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)
        self.fields = fields
        self.out_checkpoint_file = XDMFCheckpointHandler(file_name=out_checkpoint_file_name, append_to_existing=False,
                                                         field=self.fields.u_new,
                                                         field_name=self.fields.u_new.name())
        self.dddb_loader = DDDbLoader(database_file_name=dddb_file_name, strains_file_prefix='strain', params_file_prefix='params',
                                      iteration_number=10, path_to_files=None)
        self.params_db = self.dddb_loader.parameters_array
        self.strains_db = self.dddb_loader.strains_array
        self.parameter_updates = ParameterUpdates(fields=self.fields, strains_db=self.strains_db,
                                                 params_db=self.params_db, old_indices=self.random_parameter_indices())

    def random_parameter_indices(self) -> np.ndarray:
        return np.random.choice(range(len(self.params_db)), int(self.fields.function_space.dim()))

    def norm_max(self, matrix):
        return np.max(np.abs(matrix))

    def norm_sqrd(self, matrix):
        return np.sum(np.square(matrix))


    def run(self, i):
        it = 0
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t_float, i)
        while True:
            it += 1
            self.fem_solver.run(fields=self.fields)
            
            self.parameter_updates.run(norm=self.norm_sqrd)

            if sum(self.parameter_updates.new_indices == self.parameter_updates.old_indices)/len(self.parameter_updates.new_indices) > 0.9:
                print('exit')
                break
            if it % 100 == 0:
                print(sum(self.parameter_updates.new_indices == self.parameter_updates.old_indices)/len(self.parameter_updates.new_indices))
                print(all(self.parameter_updates.new_indices == self.parameter_updates.old_indices))
            self.parameter_updates.set_next_state()

        self.out_checkpoint_file.write(i)
        self.field_updates.run(fields=self.fields)

    def close(self):
        self.out_checkpoint_file.close()

