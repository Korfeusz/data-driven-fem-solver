from file_handling import HDF5File, XDMFCheckpointHandler
from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from fields.field_updates import FieldUpdates
from fields.elastodynamics_fields import ElastodynamicsFields


class ElastodynamicsTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: ElastodynamicsFields,
                 checkpoint_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)
        self.checkpoint_file = XDMFCheckpointHandler(file_name=checkpoint_file_name, append_to_existing=False,
                                                     field=fields.u_new, field_name=fields.u_new.name())



    def run(self, i: int):
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t_float, i)
        self.fem_solver.run(self.fields)
        self.field_updates.run(fields=self.fields)
        self.checkpoint_file.write(i)

    def close(self):
        self.checkpoint_file.close()
