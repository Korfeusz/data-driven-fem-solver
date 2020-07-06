from file_handling import HDF5File
from space_definition import Spaces
from .time_step import TimeStep
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from fem_solver import FemSolver
import fenics
from problem_definition.external_excitation import ExternalExcitation
from problem_definition.field_updates import FieldUpdates
from problem_definition.elastodynamics_fields import ElastodynamicsFields
import numpy as np



class DDDbTimeStep(TimeStep):
    def __init__(self, alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 fem_solver: FemSolver,
                 file: fenics.XDMFFile,
                 boundary_excitation: ExternalExcitation,
                 field_updates: FieldUpdates,
                 fields: ElastodynamicsFields,
                 mesh: fenics.Mesh,
                 hdf_file_name: str,
                 spaces: Spaces,
                 strain_file_name: str,
                 material_parameters_file_name: str):
        super().__init__(alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields)
        self.hdf5file = HDF5File(mesh=mesh, mode='r', file_name=hdf_file_name,
                                 function=fields.u_new, function_name=fields.u_new.name())
        self.tensor_space = spaces.tensor_space
        # self.t2d = fenics.vertex_to_dof_map(spaces.tensor_space)
        self.number_of_vertices = mesh.num_vertices()


    def transform_material_parameters(self, u_new: fenics.Function,
                                      new_constitutive_relation_multiplicative_parameters: fenics.Function):
        strain_vec = fenics.project(fenics.sym(fenics.grad(u_new)), V=self.tensor_space).vector()[:]
        strain_tens = np.reshape(strain_vec, newshape=(self.number_of_vertices, 2, 2))
        parameters_vec = fenics.project(new_constitutive_relation_multiplicative_parameters,
                                        V=self.tensor_space).vector()[:]
        parameters_tens = np.reshape(parameters_vec, newshape=(self.number_of_vertices, 2, 2))
        return strain_tens, parameters_tens


    def run(self, i: int):
        pass

    def close(self):
        self.hdf5file.close()
