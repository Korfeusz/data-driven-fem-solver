from typing import Type

from constitutive_relations import LinearHookesLaw
from fem_solver import FemSolver, PreAssembledSolver
from problem_definition import TimeStep, Fields, ElastodynamicsFields, ElastodynamicsTimeStep
from space_definition import Spaces, VectorFunctionSpaceCreator, TensorFunctionSpaceCreator
from .simulation_parameters import SimulationParameters


class ElastodynamicSimulationParameters(SimulationParameters):
    def __init__(self):
        self._constitutive_relation = None


    @property
    def spaces(self) -> Spaces:
        vector_space_creator = VectorFunctionSpaceCreator(element_family='CG', degree=1)
        tensor_space_creator = TensorFunctionSpaceCreator(element_family='CG', degree=1)
        return Spaces(vector_space_creator=vector_space_creator, tensor_space_creator=tensor_space_creator)


    @property
    def fields(self) -> Fields:
        return ElastodynamicsFields()


    @property
    def fem_solver_type(self) -> Type[FemSolver]:
        return PreAssembledSolver

    @property
    def constitutive_relation(self):
        if self._constitutive_relation is None:
            self._constitutive_relation = LinearHookesLaw(young_modulus=1000.0, poisson_coefficient=0.3)
        return self._constitutive_relation


    @property
    def time_step_type(self) -> Type[TimeStep]:
        return ElastodynamicsTimeStep


    @property
    def save_file_name(self) -> str:
        return 'test_results.xdmf'

