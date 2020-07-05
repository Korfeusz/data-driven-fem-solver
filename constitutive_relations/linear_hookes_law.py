from .constitutive_relation import ConstitutiveRelation
import fenics


class LinearHookesLaw(ConstitutiveRelation):
    def __init__(self, young_modulus: float, poisson_coefficient: float):
        self.young_modulus = young_modulus
        self.poisson_coefficient = poisson_coefficient

    @property
    def lame_mu(self) -> fenics.Constant:
        return fenics.Constant(self.young_modulus / (2.0 * (1.0 + self.poisson_coefficient)))

    @property
    def lame_lambda(self) -> fenics.Constant:
        return fenics.Constant(self.young_modulus*self.poisson_coefficient / ((1.0 + self.poisson_coefficient)
                                                                              *(1.0 - 2.0*self.poisson_coefficient)))

    def get_new_value(self, r: fenics.Function) -> fenics.Expression:
        return 2.0*self.lame_mu*fenics.sym(fenics.grad(r)) \
               + self.lame_lambda*fenics.tr(fenics.sym(fenics.grad(r)))*fenics.Identity(len(r))

    def get_old_value(self, r: fenics.Function) -> fenics.Expression:
        return self.get_new_value(r)