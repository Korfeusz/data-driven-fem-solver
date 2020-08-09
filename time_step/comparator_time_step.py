from typing import Optional

from file_handling import HDF5File, XDMFCheckpointHandler, NPYFile
from optimizer import ScipyOptimizer
from fields import DDDbFields
from space_definition import Spaces
from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from fields.field_updates import FieldUpdates
import numpy as np



class ComparatorTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: DDDbFields,
                 mesh: fenics.Mesh,
                 in_file_iteration_number: Optional[int],
                 in_checkpoint_file_name: str,
                 out_checkpoint_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)


    def run(self, i):
        pass

    def close(self):
        pass
