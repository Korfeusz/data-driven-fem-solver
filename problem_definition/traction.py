from .external_excitation import ExternalExcitation
import fenics
from generalized_alpha_parameters import GeneralizedAlphaParameters


class Traction(ExternalExcitation):
    def __init__(self, initial_value: float, cutoff_time: float, traction_boundary: str):
        self.force = fenics.Expression(("t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_time,
                                       p0=initial_value, degree=1)
        self.traction_boundary = traction_boundary
        self.ds = None

    def set_ds(self, ds: fenics.Measure):
        self.ds = ds

    def value(self, w: fenics.Function):
        return fenics.dot(w, self.force)*self.ds(self.traction_boundary)

    def update(self, alpha_params: GeneralizedAlphaParameters, delta_t: float, iteration: int):
        # self.force.t = (iteration + 1. - alpha_params.alpha_f_float) * delta_t
        self.force.t = iteration * delta_t