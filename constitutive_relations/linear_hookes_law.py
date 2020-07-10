from typing import Union

from lazy import Lazy
from fields import DDDbFields
from .constitutive_relation import ConstitutiveRelation
import fenics


class LinearHookesLaw(ConstitutiveRelation):
    def __init__(self, poisson_coefficient: float, young_modulus_or_fields: Union[float, DDDbFields]):
        self._young_modulus_or_fields = young_modulus_or_fields
        self.poisson_coefficient = fenics.Constant(poisson_coefficient)

    @property
    @Lazy
    def young_modulus(self) -> fenics.Expression:
        if isinstance(self._young_modulus_or_fields, DDDbFields):
            young_modulus = self._young_modulus_or_fields.new_constitutive_relation_multiplicative_parameters
        elif isinstance(self._young_modulus_or_fields, float):
            young_modulus = fenics.Constant(self._young_modulus_or_fields)
        else:
            raise TypeError('Please specify either a float or a field for the LinearHookesLaw young modulus parameter')
        return young_modulus

    @property
    def lame_mu(self) -> fenics.Expression:
        return self.young_modulus / fenics.Constant(2.0 * (1.0 + self.poisson_coefficient))

    @property
    def lame_lambda(self) -> fenics.Expression:
        return self.young_modulus*fenics.Constant(self.poisson_coefficient / ((1.0 + self.poisson_coefficient)
                                                                              *(1.0 - 2.0*self.poisson_coefficient)))

    def get_new_value(self, r: fenics.Function) -> fenics.Expression:
        return 2.0*self.lame_mu*fenics.sym(fenics.grad(r)) \
               + self.lame_lambda*fenics.tr(fenics.sym(fenics.grad(r)))*fenics.Identity(len(r))

    def get_old_value(self, r: fenics.Function) -> fenics.Expression:
        return self.get_new_value(r)