import abc
from typing import Type

class FemSolver(abc.ABC):
    @abc.abstractmethod
    def __init__(self, problem, fields, boundary_conditions):
        pass

    @abc.abstractmethod
    def run(self, fields):
        pass


def get_fem_solver(fem_solver: Type[FemSolver], problem, fields, boundary_conditions):
    return fem_solver(problem, fields, boundary_conditions)