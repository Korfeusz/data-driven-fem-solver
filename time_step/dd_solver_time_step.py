import fenics

from fem_solver import FemSolver
from fields import FieldUpdates, DDDbFields
from file_handling import XDMFCheckpointHandler
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation
from time_step import TimeStep
from time_stepping_parameters import TimeSteppingParameters
from data_driven_solver import DDDbLoader
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

    def random_parameter_indices(self) -> np.ndarray:
        return np.random.choice(range(len(self.params_db)), int(self.fields.function_space.dim()))

    def create_strain_tensor(self, u):
        strain_vec = fenics.project(fenics.sym(fenics.grad((u))), V=self.fields.tensor_space).vector()[:]
        return np.reshape(strain_vec, newshape=(self.fields.function_space.dim(), 2, 2))

    def norm_max(self, matrix):
        return np.max(np.abs(matrix))

    def norm_sqrd(self, matrix):
        return np.sum(np.square(matrix))

    def get_index_of_closest_matrix(self, matrix, matrix_set, norm):
        distances = [norm(x - matrix) for x in matrix_set]
        return np.argmin(distances)

    def get_new_young_moduli_indices(self, computed_strains, strain_dataset, norm):
        indices = np.zeros(len(computed_strains), dtype=int)
        for i, strain in enumerate(computed_strains):
            closest_index = self.get_index_of_closest_matrix(strain, strain_dataset, norm)
            indices[i] = closest_index
        return indices

    def run(self, i):
        old_indices = self.random_parameter_indices()
        parameters = self.params_db[old_indices]
        strains = self.strains_db[old_indices]
        self.fields.new_constitutive_relation_multiplicative_parameters = parameters
        it = 0
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t_float, i)
        while True:
            it += 1
            self.fem_solver.run(fields=self.fields)

            strain_tens = self.create_strain_tensor(self.fields.u_new)
            new_indices = self.get_new_young_moduli_indices(strain_tens, self.strains_db, self.norm_sqrd)
            self.fields.new_constitutive_relation_multiplicative_parameters = self.params_db[new_indices]

            # new_strains = np.reshape(self.strains_db[new_indices], newshape=(self.fields.tensor_space.dim()))
            # strain_tensor.vector()[:] = new_strains
            if sum(new_indices == old_indices)/len(new_indices) > 0.9:
                print('exit')
                break
            if it % 100 == 0:
                print(sum(new_indices == old_indices)/len(new_indices))
                print(all(new_indices == old_indices))
            old_indices = new_indices

        self.out_checkpoint_file.write(i)
        self.field_updates.run(fields=self.fields)

    def close(self):
        self.out_checkpoint_file.close()

