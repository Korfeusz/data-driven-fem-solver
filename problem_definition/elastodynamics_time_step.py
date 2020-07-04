from .time_step import TimeStep


class ElastodynamicsTimeStep(TimeStep):
    def __init__(self, alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields):
        super().__init__(alpha_params, time_params, fem_solver, file, boundary_excitation, field_updates, fields)


    def run(self, i):
        self.boundary_excitation.update(self.alpha_params, self.time_params.delta_t, i)
        self.fem_solver.run(self.fields)
        self.field_updates.run(fields=self.fields)
        self.file.write(self.fields.u_new, (i + 1)*self.time_params.delta_t_float)
