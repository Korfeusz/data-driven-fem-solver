from typing import Type

from fem_solver import FemSolver
from fields import Fields
from space_definition import Spaces
from time_step import TimeStepBuilder
from .simulation_parameters import SimulationParameters

class DDSolverSimulationParameters(SimulationParameters):
    def __init__(self):
        pass

    @property
    def spaces(self) -> Spaces:
        pass

    @property
    def fields(self) -> Fields:
        pass

    @property
    def fem_solver_type(self) -> Type[FemSolver]:
        pass

    @property
    def time_step_builder(self) -> TimeStepBuilder:
        pass

    @property
    def constitutive_relation(self):
        pass