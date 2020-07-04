from space_definition import  BoundaryMarkers, DirichletBCCreator, Spaces, MeshCreator
import fenics
from problem_definition import  TimeStepBuilder, ExternalExcitation, Fields, ProblemForm, TimeStep, FieldUpdates
from fem_solver import get_fem_solver, FemSolver
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters
from typing import Type


class Simulation:
    def __init__(self,
                 mesh_creator: MeshCreator,
                 spaces: Spaces,
                 boundary_markers: BoundaryMarkers,
                 bc_creator: DirichletBCCreator,
                 boundary_excitation: ExternalExcitation,
                 fields: Fields,
                 field_updates: FieldUpdates,
                 fem_solver_type: Type[FemSolver],
                 problem: ProblemForm,
                 time_step_type: Type[TimeStep],
                 alpha_params: GeneralizedAlphaParameters,
                 time_params: TimeSteppingParameters,
                 save_file_name: str
                 ):
        self.mesh_creator = mesh_creator
        self.spaces = spaces
        self.boundary_markers = boundary_markers
        self.bc_creator = bc_creator
        self.boundary_excitation = boundary_excitation
        self.fields = fields
        self.field_updates = field_updates
        self.fem_solver_type = fem_solver_type
        self.problem = problem
        self.alpha_params = alpha_params
        self.time_params = time_params
        self.time_step_builder = TimeStepBuilder(time_step_type=time_step_type)
        self.xdmf_file = fenics.XDMFFile(save_file_name)
        self.xdmf_file.parameters["flush_output"] = True
        self.xdmf_file.parameters["functions_share_mesh"] = True
        self.xdmf_file.parameters["rewrite_function_mesh"] = False

    def run(self):
        mesh = self.mesh_creator.get_mesh()
        self.spaces.generate(mesh=mesh)
        self.boundary_markers.mark_boundaries(mesh=mesh)
        bc = self.bc_creator.apply(vector_space=self.spaces.vector_space, boundary_markers=self.boundary_markers.value)
        ds = fenics.Measure('ds', domain=mesh, subdomain_data=self.boundary_markers.value)
        self.boundary_excitation.set_ds(ds=ds)
        self.fields.generate(spaces=self.spaces)
        fem_solver = get_fem_solver(fem_solver=self.fem_solver_type, problem=self.problem, fields=self.fields,
                                    boundary_conditions=bc)
        self.time_step_builder.set(alpha_params=self.alpha_params, time_params=self.time_params, fem_solver=fem_solver,
                              file=self.xdmf_file, boundary_excitation=self.boundary_excitation,
                              field_updates=self.field_updates, fields=self.fields)
        time_step = self.time_step_builder.build()

        for (i, t) in enumerate(self.time_params.linear_time_space[1:]):
            print("Time: ", t)
            time_step.run(i)