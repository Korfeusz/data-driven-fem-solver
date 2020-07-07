import fenics

from problem_definition import DDDbFields
from .optimizer import Optimizer
import scipy.optimize as spo
import numpy as np
from fem_solver import FemSolver

class ScipyOptimizer(Optimizer):
    def __init__(self, fields: DDDbFields, initial_values: np.ndarray, fem_solver: FemSolver):
        self.parameters = initial_values
        self.fields = fields
        self.field_to_optimize = fields.new_constitutive_relation_multiplicative_parameters
        self.fem_solver = fem_solver

    def compare_displacement(self, u_1, u_2, vector_space, function_space):
        diff = fenics.project(u_1 - u_2, V=vector_space)
        u_norm = np.zeros(fs.dim())
        for f in fenics.split(diff):
            u_norm += fenics.project(f, V=function_space).vector()[v2d] ** 2
        u_norm = np.sqrt(u_norm)
        return np.sqrt(np.sum(np.power(u_norm, 2)))

    def function_to_minimize(self, parameters):
        self.field_to_optimize.vector()[:] = parameters
        self.fem_solver.run(self.fields)
        return self.compare_displacement()

    def run(self):
        self._results = spo.minimize(self.function_to_minimize, self.parameters, args=(a, L, v2d, bc, u, u_old, time_generator()),
                                      tol=1e-6, method='SLSQP', bounds=bounds, options={'maxiter': 1e6})

    @property
    def results(self):
        return self._results