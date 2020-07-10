import abc
from typing import Type

from fem_solver import FemSolver
from fields import Fields
from time_step import TimeStepBuilder
from space_definition import Spaces


class SimulationParameters(abc.ABC):
    @property
    @abc.abstractmethod
    def spaces(self) -> Spaces:
        pass

    @property
    @abc.abstractmethod
    def fields(self) -> Fields:
        pass


    @property
    @abc.abstractmethod
    def fem_solver_type(self) -> Type[FemSolver]:
        pass


    @property
    @abc.abstractmethod
    def time_step_builder(self) -> TimeStepBuilder:
        pass


    @property
    @abc.abstractmethod
    def save_file_name(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def constitutive_relation(self):
        pass

