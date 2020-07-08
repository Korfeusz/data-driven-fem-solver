import fenics
from problem_definition import DDDbFields
from .optimizer import Optimizer
import scipy.optimize as spo
import numpy as np
from fem_solver import FemSolver
from space_definition import Spaces

class ScipyOptimizer(Optimizer):
    def __init__(self, fields: DDDbFields, initial_values: np.ndarray, fem_solver: FemSolver, spaces: Spaces):
        self.parameters = (initial_values * 10) - 5.
        self.fields = fields
        # self.field_to_optimize = fields.new_constitutive_relation_multiplicative_parameters
        self.fem_solver = fem_solver
        self.spaces = spaces
        self._results = None
        self.it = 0

    def compare_displacement(self, fields):
        diff = fenics.project(fields.u_new - fields.imported_displacement_field, V=self.spaces.vector_space)
        u_norm = np.zeros(self.spaces.function_space.dim())
        for f in fenics.split(diff):
            u_norm += np.power(fenics.project(f, V=self.spaces.function_space).vector()[:], 2)
        u_norm = np.sqrt(u_norm)
        norm = np.sqrt(np.sum(np.power(u_norm, 2)))
        # print('norm: {}'.format(norm))
        return norm

    def function_to_minimize(self, parameters):
        self.fields.new_constitutive_relation_multiplicative_parameters.vector()[:] = parameters
        self.fem_solver.run(self.fields)
        result = self.compare_displacement(self.fields)
        if self.it % 100 == 0:
            print('iteration: {}, norm: {}'.format(self.it, result))
        self.it +=1
        return result

    def run(self):
        # for i in range(3):
        #     self.fields.new_constitutive_relation_multiplicative_parameters.vector()[:] = np.random.random(self.spaces.function_space.dim() * 4)*10
        #     self.fem_solver.run(self.fields)
        #     print(self.fields.u_new.vector()[0:10])
        #     self.compare_displacement(self.fields)
        self.it = 0
        # bounds = tuple((0., 100.) for _ in range(self.spaces.function_space.dim() * 4))
        bounds = None
        self._results = spo.minimize(self.function_to_minimize, self.parameters, tol=1e-12, method='SLSQP', bounds=bounds,
                                     options={'maxiter': 1e6})

    @property
    def results(self):
        return self._results