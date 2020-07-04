import abc
from space_definition import  BoundaryMarkers, DirichletBCCreator, Spaces, MeshCreator
from problem_definition import  ExternalExcitation, Fields, ProblemForm, TimeStep, FieldUpdates
from fem_solver import FemSolver
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from typing import Type



class SimulationParameters(abc.ABC):
    @property
    @abc.abstractmethod
    def mesh_creator(self) -> MeshCreator:
        pass

    @property
    @abc.abstractmethod
    def spaces(self) -> Spaces:
        pass

    @property
    @abc.abstractmethod
    def boundary_markers(self) -> BoundaryMarkers:
        pass

    @property
    @abc.abstractmethod
    def bc_creator(self) -> DirichletBCCreator:
        pass

    @property
    @abc.abstractmethod
    def boundary_excitation(self) -> ExternalExcitation:
        pass

    @property
    @abc.abstractmethod
    def fields(self) -> Fields:
        pass

    @property
    @abc.abstractmethod
    def field_updates(self) -> FieldUpdates:
        pass

    @property
    @abc.abstractmethod
    def fem_solver_type(self) -> Type[FemSolver]:
        pass

    @property
    @abc.abstractmethod
    def problem(self) -> ProblemForm:
        pass

    @property
    @abc.abstractmethod
    def time_step_type(self) -> Type[TimeStep]:
        pass

    @property
    @abc.abstractmethod
    def alpha_params(self) -> GeneralizedAlphaParameters:
        pass

    @property
    @abc.abstractmethod
    def time_params(self) -> TimeSteppingParameters:
        pass

    @property
    @abc.abstractmethod
    def save_file_name(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def constitutive_relation(self):
        pass

