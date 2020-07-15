from fem_solver import FemSolver
from fields import FieldUpdates, DDDbFields
from file_handling import XDMFCheckpointHandler
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation
from time_step import TimeStep
from time_stepping_parameters import TimeSteppingParameters
from data_driven_solver import DDDbLoader


class DDSolverTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: DDDbFields,
                 out_checkpoint_file_name: str,
                 dddb_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, boundary_excitation, field_updates, fields)
        self.out_checkpoint_file = XDMFCheckpointHandler(file_name=out_checkpoint_file_name, append_to_existing=False,
                                                         field=fields.u_new,
                                                         field_name=fields.u_new.name())
        self.dddb_loader = DDDbLoader(database_file_name=dddb_file_name, strains_file_prefix='strains', params_file_prefix='params',
                                      iteration_number=10, path_to_files=None)
        self.params_db = self.dddb_loader.parameters_array
        self.strains_db = self.dddb_loader.strains_array


    def random_parameter_initialization(self):
        self.dddb_loader.parameters_array


    def run(self, i):
        pass

    def close(self):
        pass

