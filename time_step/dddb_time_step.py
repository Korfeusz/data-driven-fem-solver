from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from problem_definition.field_updates import FieldUpdates
from problem_definition.elastodynamics_fields import ElastodynamicsFields

class DDDbTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 file: fenics.XDMFFile,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: ElastodynamicsFields):
        super().__init__(alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields)


    def run(self, i: int):
        pass
