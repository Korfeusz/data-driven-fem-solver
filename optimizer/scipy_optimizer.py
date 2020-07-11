import fenics
from fields import DDDbFields
from .optimizer import Optimizer
import scipy.optimize as spo
import numpy as np
from fem_solver import FemSolver
from space_definition import Spaces

class ScipyOptimizer(Optimizer):
    def __init__(self, fields: DDDbFields, initial_values: np.ndarray, fem_solver: FemSolver, spaces: Spaces):
        self.parameters = (initial_values * 10) - 5.
        self.fields = fields
        self.fem_solver = fem_solver
        self.spaces = spaces
        self._results = None
        self.it = 0
        self.exit_value = 1e-2
        self.norm = None

    def exit_condition(self, xk):
        if self.norm < self.exit_value:
            raise RuntimeError

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
        self.norm = result
        if self.it % 100 == 0:
            print('iteration: {}, norm: {}'.format(self.it, result))
        self.it +=1
        return result


    def run(self):
        self.it = 0
        bounds = None
        try:
            self._results = spo.minimize(self.function_to_minimize, self.parameters, tol=1e-15, method='SLSQP', bounds=bounds,
                                     options={'maxiter': 1e6}, callback=self.exit_condition)
        except RuntimeError:
            print('Exited optimizer, condition met')
            pass
    #         save parameters and fields here

    @property
    def results(self):
        return self._results