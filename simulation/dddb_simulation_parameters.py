from typing import Type

from constitutive_relations import ConstitutiveRelation
from constitutive_relations.dddb_constitutive_relation import DDDbConstitutiveRelation
from fem_solver import FemSolver, PreAssembledSolver
from problem_definition import Fields
from problem_definition.dddb_fields import DDDbFields
from space_definition import Spaces, VectorFunctionSpaceCreator, TensorFunctionSpaceCreator
from time_step import TimeStep, TimeStepBuilder
from time_step.dddb_time_step import DDDbTimeStep
from .simulation_parameters import SimulationParameters
import numpy as np

class DDDbParameters(SimulationParameters):
    def __init__(self):
        self._constitutive_relation = None
        self._fields = None


    @property
    def spaces(self) -> Spaces:
        vector_space_creator = VectorFunctionSpaceCreator(element_family='CG', degree=1)
        tensor_space_creator = TensorFunctionSpaceCreator(element_family='CG', degree=1)
        return Spaces(vector_space_creator=vector_space_creator, tensor_space_creator=tensor_space_creator)


    @property
    def fields(self) -> DDDbFields:
        if self._fields is None:
            self._fields = DDDbFields()
        return self._fields


    @property
    def fem_solver_type(self) -> Type[FemSolver]:
        return PreAssembledSolver

    @property
    def constitutive_relation(self) -> ConstitutiveRelation:
        if self._constitutive_relation is None:
            self._constitutive_relation = DDDbConstitutiveRelation(fields=self.fields)
        return self._constitutive_relation


    @property
    def time_step_builder(self) -> TimeStepBuilder:
        tsb = TimeStepBuilder(time_step_type=DDDbTimeStep)
        tsb.set(hdf5_file_name='saving_elastic.h5',
                material_parameters_file_name='material_params.npy',
                strain_file_name='strain.npy')
        return tsb


    @property
    def save_file_name(self) -> str:
        return 'db_test_results.xdmf'

