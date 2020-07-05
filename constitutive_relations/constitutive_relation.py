import abc
import fenics

class ConstitutiveRelation(abc.ABC):
    @abc.abstractmethod
    def get_new_value(self, r: fenics.Function) -> fenics.Expression:
        pass

    @abc.abstractmethod
    def get_old_value(self, r: fenics.Function) -> fenics.Expression:
        pass