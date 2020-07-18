import numpy as np
from .norm_builder import NormBuilder
import fenics


class ParameterUpdates:
    def __init__(self, fields, strains_db, params_db, old_indices):
        self.fields = fields
        self.strains_db = strains_db
        self.params_db = params_db
        self.new_indices = None
        self.old_indices = old_indices
        # params transformation to the correct shape
        self.fields.new_constitutive_relation_multiplicative_parameters = self.params_db[old_indices]


    def run(self, norm):
        computed_strains = self.create_strain_tensor(self.fields.u_new)
        self.new_indices = self.get_new_indices(computed_strains, norm)
        self.fields.new_constitutive_relation_multiplicative_parameters = self.params_db[self.new_indices]

    def set_next_state(self):
        self.old_indices = self.new_indices

    def get_index_of_closest_matrix(self, matrix, norm):
        distances = [norm(x - matrix) for x in self.strains_db]
        return np.argmin(distances)

    def get_new_indices(self, computed_strains, norm):
        indices = np.zeros(len(computed_strains), dtype=int)
        for i, strain in enumerate(computed_strains):
            closest_index = self.get_index_of_closest_matrix(strain, norm)
            indices[i] = closest_index
        return indices

    def create_strain_tensor(self, u):
        strain_vec = fenics.project(fenics.sym(fenics.grad((u))), V=self.fields.tensor_space).vector()[:]
        return np.reshape(strain_vec, newshape=(self.fields.function_space.dim(), 2, 2))


