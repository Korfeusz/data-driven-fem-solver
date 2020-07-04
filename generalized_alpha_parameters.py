import fenics
from typing import Optional

class GeneralizedAlphaParameters:
    def __init__(self, alpha_m: float, alpha_f:float, gamma: Optional[float] = None, beta: Optional[float] = None):
        self.alpha_m = fenics.Constant(alpha_m)
        self.alpha_f = fenics.Constant(alpha_f)
        self.alpha_f_float = alpha_f
        if gamma is None:
            self.gamma = fenics.Constant(0.5+alpha_f-alpha_m)
        else:
            self.gamma = fenics.Constant(gamma)
        if beta is None:
            self.beta = fenics.Constant((self.gamma+0.5)**2/4.)
        else:
            self.beta = fenics.Constant(beta)