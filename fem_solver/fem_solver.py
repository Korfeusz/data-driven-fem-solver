import abc
from typing import Type, List
from problem_definition import ProblemForm
from fields import Fields
import fenics

class FemSolver(abc.ABC):
    @abc.abstractmethod
    def __init__(self, problem: ProblemForm, fields: Fields, boundary_conditions: List[fenics.DirichletBC]):
        pass

    @abc.abstractmethod
    def run(self, fields: Fields):
        pass


def get_fem_solver(fem_solver: Type[FemSolver], problem: ProblemForm, fields: Fields, boundary_conditions: List[fenics.DirichletBC]):
    return fem_solver(problem, fields, boundary_conditions)