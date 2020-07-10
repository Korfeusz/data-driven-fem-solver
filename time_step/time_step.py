import abc
import fenics
from fem_solver import FemSolver
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation
from fields import FieldUpdates
from time_stepping_parameters import TimeSteppingParameters


class TimeStep(abc.ABC):
    @abc.abstractmethod
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 file: fenics.XDMFFile,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields):
        self.alpha_params = alpha_params
        self.time_params = time_params
        self.fem_solver = fem_solver
        self.file = file
        self.boundary_excitation = boundary_excitation
        self.field_updates = field_updates
        self.fields = fields

    @abc.abstractmethod
    def run(self, i):
        pass

    @abc.abstractmethod
    def close(self):
        pass