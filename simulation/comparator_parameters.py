from typing import Type, Optional
from constitutive_relations import ConstitutiveRelation
from fem_solver import FemSolver, PreAssembledSolver
from fields import Fields
from lazy import Lazy
from space_definition import Spaces, VectorFunctionSpaceCreator, TensorFunctionSpaceCreator
from space_definition.normal_function_space_creator import NormalFunctionSpaceCreator
from time_step import TimeStepBuilder, ComparatorTimeStep
from .simulation_parameters import SimulationParameters


class ComparatorParameters(SimulationParameters):
    def __init__(self):
       pass


    @property
    def spaces(self) -> Spaces:
        vector_space_creator = VectorFunctionSpaceCreator(element_family='CG', degree=1)
        tensor_space_creator = TensorFunctionSpaceCreator(element_family='CG', degree=1)
        function_space_creator = NormalFunctionSpaceCreator(element_family='CG', degree=1)
        return Spaces(vector_space_creator=vector_space_creator, tensor_space_creator=tensor_space_creator,
                      function_space_creator=function_space_creator)


    @property
    @Lazy
    def fields(self) -> Fields:
        return Fields()


    @property
    def fem_solver_type(self) -> Type[FemSolver]:
        return PreAssembledSolver

    @property
    @Lazy
    def constitutive_relation(self) -> Optional[ConstitutiveRelation]:
        return None


    @property
    def time_step_builder(self) -> TimeStepBuilder:
        tsb = TimeStepBuilder(time_step_type=ComparatorTimeStep)
        tsb.set(in_checkpoint_file_name='checkpoint_file.xdmf',
                out_checkpoint_file_name='dddb_out_file.xdmf',
                dddb_file_name='dddb_out_tensor',
                )
        return tsb
