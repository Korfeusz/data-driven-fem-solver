from file_handling import HDF5File
from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from problem_definition.field_updates import FieldUpdates
from problem_definition.elastodynamics_fields import ElastodynamicsFields

class ElastodynamicsTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 file: fenics.XDMFFile,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: ElastodynamicsFields,
                 mesh: fenics.Mesh,
                 hdf_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields)

        self.hdf5file = HDF5File(mesh=mesh, mode='w',file_name=hdf_file_name,
                                 function=fields.u_new, function_name=fields.u_new.name())



    def run(self, i: int):
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t, i)
        self.fem_solver.run(self.fields)
        self.field_updates.run(fields=self.fields)
        self.file.write(self.fields.u_new, (i + 1)*self.time_params.delta_t_float)
        self.hdf5file.save(i)
