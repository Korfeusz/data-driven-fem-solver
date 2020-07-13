from typing import Optional

import fenics
from fields import DDDbFields
from .optimizer import Optimizer
import scipy.optimize as spo
import numpy as np
from fem_solver import FemSolver
from space_definition import Spaces
import threading
from file_handling import XDMFCheckpointHandler

class ScipyOptimizer(Optimizer):
    def __init__(self, fields: DDDbFields, initial_values: np.ndarray, fem_solver: FemSolver, spaces: Spaces):
        self.parameters = (initial_values * 20.) + 990.
        self.fields = fields
        self.fem_solver = fem_solver
        self.spaces = spaces
        self._results = None
        self.it = 0
        self.exit_value = 1e-6
        self.norm = None
        self.keep_going = True
        self.stop_thread = threading.Thread(target=self.key_capture_thread, args=(), name='key_capture_thread', daemon=True).start()
        self.save_file: Optional[XDMFCheckpointHandler] = None

    def key_capture_thread(self):
        input()
        self.keep_going = False

    def optimizer_callback(self, xk):
        if self.it % 1 == 0:
            print('iteration: {}, norm: {}'.format(self.it, self.norm))
            if self.save_file is not None:
                self.save_file.write(self.it)
        self.it +=1
        if self.norm < self.exit_value:
            print('Condition met')
            raise RuntimeError
        elif not self.keep_going:
            print('User exit')
            raise RuntimeError
        # If it does not change enough between steps shake it up a little

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
        self.norm = self.compare_displacement(self.fields)
        return self.norm


    def run(self):
        self.it = 0
        bounds = tuple((990., 1100.) for _ in range(self.fields.new_constitutive_relation_multiplicative_parameters.function_space().dim()))

        try:
            spo.minimize(self.function_to_minimize, self.parameters, tol=0.0, method='SLSQP', bounds=bounds,
                         options={'maxiter': 1e6, 'ftol': 0.0}, callback=self.optimizer_callback)
            # spo.minimize(self.function_to_minimize, self.parameters, tol=0.0, method='Powell',
            #              options={'maxiter': 1e6, 'xtol': 0.0, 'ftol': 0.0, 'disp': True}, callback=self.optimizer_callback)
        except RuntimeError:
            print('Exited optimizer')
            pass
