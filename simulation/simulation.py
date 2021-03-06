import fenics
from fem_solver import get_fem_solver
from .simulation_parameters import SimulationParameters
from .common_simulation_parameters import CommonSimulationParameters
import numpy as np

class Simulation:
    def __init__(self, simulation_parameters: SimulationParameters,
                 common_simulation_parameters: CommonSimulationParameters):
        self.mesh_creator = common_simulation_parameters.mesh_creator
        self.spaces = simulation_parameters.spaces
        self.boundary_markers = common_simulation_parameters.boundary_markers
        self.bc_creator = common_simulation_parameters.bc_creator
        self.boundary_excitation = common_simulation_parameters.boundary_excitation
        self.fields = simulation_parameters.fields
        self.field_updates = common_simulation_parameters.field_updates
        self.fem_solver_type = simulation_parameters.fem_solver_type
        self.problem = common_simulation_parameters.problem(simulation_parameters.constitutive_relation)
        self.alpha_params = common_simulation_parameters.alpha_params
        self.time_params = common_simulation_parameters.time_params
        self.time_step_builder = simulation_parameters.time_step_builder



    def run(self) -> None:
        mesh = self.mesh_creator.get_mesh()
        self.spaces.initialize(mesh=mesh)
        self.boundary_markers.mark_boundaries(mesh=mesh)
        bc = self.bc_creator.apply(vector_space=self.spaces.vector_space, boundary_markers=self.boundary_markers.value)
        ds = fenics.Measure('ds', domain=mesh, subdomain_data=self.boundary_markers.value)
        self.boundary_excitation.set_ds(ds=ds)
        self.fields.initialize(spaces=self.spaces)
        fem_solver = get_fem_solver(fem_solver=self.fem_solver_type, problem=self.problem, fields=self.fields,
                                    boundary_conditions=bc)
        self.time_step_builder.set(alpha_params=self.alpha_params, time_params=self.time_params, fem_solver=fem_solver,
                              boundary_excitation=self.boundary_excitation,
                              field_updates=self.field_updates, fields=self.fields, mesh=mesh, spaces=self.spaces)
        time_step = self.time_step_builder.build()
        mins = []
        maxs = []
        for (i, t) in enumerate(self.time_params.linear_time_space[1:]):
        # for i in range(3):
            print("Time: ", t)
            time_step.run(i)
            # v = fenics.elem_div(self.problem.constitutive_relation.get_new_value(self.fields.u_new), fenics.sym(fenics.grad(self.fields.u_new)))
            # x = fenics.project(v, self.spaces.tensor_space).vector()[:]
            # mins.append(min(x))
            # maxs.append(max(x))
            # # print('min: {}, max: {}, median, {}'.format(min(x), max(x), np.nanmedian(np.asarray(x))))
            # print(x)
            # print()
            # print(fenics.project(self.problem.constitutive_relation.get_new_value(self.fields.u_new), self.spaces.tensor_space).vector()[:])
            # print(fenics.project(fenics.sym(fenics.grad(self.fields.u_new)), self.spaces.tensor_space).vector()[:])
            if time_step.halt:
                break
        time_step.close()
        print('total: min: {}, max: {}'.format(min(mins), max(maxs)))