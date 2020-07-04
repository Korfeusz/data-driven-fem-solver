from .field_updates import FieldUpdates
from generalized_alpha_parameters import GeneralizedAlphaParameters
from time_stepping_parameters import TimeSteppingParameters


class ElastodynamicsFieldUpdates(FieldUpdates):
    def __init__(self, generalized_alpha_parameters: GeneralizedAlphaParameters,
                 time_stepping_parameters: TimeSteppingParameters):
        self.gamma = float(generalized_alpha_parameters.gamma)
        self.beta = float(generalized_alpha_parameters.beta)
        self.dt = time_stepping_parameters.delta_t_float

    def update_a(self, u, u_old, v_old, a_old):
        return (u - u_old - self.dt * v_old) / self.beta / self.dt ** 2 - (1 - 2 * self.beta) / 2 / self.beta * a_old

    def update_v(self, a, v_old, a_old):
        return v_old + self.dt * ((1 - self.gamma) * a_old + self.gamma * a)

    def run(self, fields):
        u = fields.u_new
        u_old = fields.u_old
        v_old = fields.v_old
        a_old = fields.a_old

        u_vec, u0_vec = u.vector(), u_old.vector()
        v0_vec, a0_vec = v_old.vector(), a_old.vector()

        # use update functions using vector arguments
        a_vec = self.update_a(u_vec, u0_vec, v0_vec, a0_vec)
        v_vec = self.update_v(a_vec, v0_vec, a0_vec)

        # Update (u_old <- u)
        v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
        u_old.vector()[:] = u.vector()
