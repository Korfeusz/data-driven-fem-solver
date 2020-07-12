from fem_solver import FemSolver
from fields import FieldUpdates
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation
from time_step import TimeStep
from time_stepping_parameters import TimeSteppingParameters


class DDSolverTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters, time_params: TimeSteppingParameters,
                 fem_solver: FemSolver, boundary_excitation: ExternalExcitation, field_updates: FieldUpdates, fields):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)

    def run(self, i):
        pass

    def close(self):
        pass

