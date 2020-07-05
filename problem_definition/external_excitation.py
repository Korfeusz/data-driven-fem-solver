import abc

import fenics

from generalized_alpha_parameters import GeneralizedAlphaParameters


class ExternalExcitation(abc.ABC):
    @abc.abstractmethod
    def value(self, w: fenics.Function) -> fenics.Expression:
        pass

    @abc.abstractmethod
    def update(self, alpha_params: GeneralizedAlphaParameters, delta_t: float, iteration: int):
        pass

    @abc.abstractmethod
    def set_ds(self, ds: fenics.Measure):
        pass
