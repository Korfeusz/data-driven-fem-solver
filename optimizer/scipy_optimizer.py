import fenics

from problem_definition import DDDbFields
from .optimizer import Optimizer
import scipy.optimize as spo
import numpy as np
from fem_solver import FemSolver
from space_definition import Spaces

class ScipyOptimizer(Optimizer):
    def __init__(self, fields: DDDbFields, initial_values: np.ndarray, fem_solver: FemSolver, spaces: Spaces):
        self.parameters = initial_values
        self.fields = fields
        self.field_to_optimize = fields.new_constitutive_relation_multiplicative_parameters
        self.fem_solver = fem_solver
        self.spaces = spaces
        self._results = None

    def compare_displacement(self, fields):
        diff = fenics.project(fields.u_new - fields.imported_displacement_field, V=self.spaces.vector_space)
        u_norm = np.zeros(self.spaces.function_space.dim())
        for f in fenics.split(diff):
            u_norm += np.power(fenics.project(f, V=self.spaces.function_space).vector()[:], 2)
        u_norm = np.sqrt(u_norm)
        return np.sqrt(np.sum(np.power(u_norm, 2)))

    def function_to_minimize(self, parameters):
        self.field_to_optimize.vector()[:] = parameters
        self.fem_solver.run(self.fields)
        return self.compare_displacement(self.fields)

    def run(self):
        bounds = None
        self._results = spo.minimize(self.function_to_minimize, self.parameters, tol=1e-6, method='SLSQP',
                                     bounds=bounds, options={'maxiter': 1e6})

    @property
    def results(self):
        return self._results