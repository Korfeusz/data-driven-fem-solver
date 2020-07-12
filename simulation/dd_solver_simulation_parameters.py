from typing import Type

from constitutive_relations import ConstitutiveRelation, LinearHookesLaw, DDDbConstitutiveRelation
from fem_solver import FemSolver, PreAssembledSolver
from fields import DDDbFields, DataDrivenParametersSpaceName
from space_definition import Spaces, VectorFunctionSpaceCreator, TensorFunctionSpaceCreator
from space_definition.normal_function_space_creator import NormalFunctionSpaceCreator
from time_step import TimeStepBuilder, DDSolverTimeStep
from .simulation_parameters import SimulationParameters

class DDSolverSimulationParameters(SimulationParameters):
    def __init__(self):
        self._fields = None
        self.constitutive_relation_type: Type[ConstitutiveRelation] = LinearHookesLaw
        self._constitutive_relation = None
        self.set_constitutive_relation_and_fields()

    def set_constitutive_relation_and_fields(self) -> None:
        if self.constitutive_relation_type is LinearHookesLaw:
            self._fields = DDDbFields(DataDrivenParametersSpaceName.FUNCTION_SPACE)
            self._constitutive_relation = LinearHookesLaw(young_modulus_or_fields=self.fields, poisson_coefficient=0.3)
        elif self.constitutive_relation_type is DDDbConstitutiveRelation:
            self._fields = DDDbFields(DataDrivenParametersSpaceName.TENSOR_SPACE)
            self._constitutive_relation = DDDbConstitutiveRelation(fields=self.fields)

    @property
    def spaces(self) -> Spaces:
        vector_space_creator = VectorFunctionSpaceCreator(element_family='CG', degree=1)
        tensor_space_creator = TensorFunctionSpaceCreator(element_family='CG', degree=1)
        function_space_creator = NormalFunctionSpaceCreator(element_family='CG', degree=1)
        return Spaces(vector_space_creator=vector_space_creator, tensor_space_creator=tensor_space_creator,
                      function_space_creator=function_space_creator)

    @property
    def fields(self) -> DDDbFields:
        return self._fields

    @property
    def fem_solver_type(self) -> Type[FemSolver]:
        return PreAssembledSolver

    @property
    def time_step_builder(self) -> TimeStepBuilder:
        tsb = TimeStepBuilder(time_step_type=DDSolverTimeStep)
        tsb.set()
        return tsb

    @property
    def constitutive_relation(self):
        return self._constitutive_relation