import abc
import fenics


class FunctionSpaceCreator(abc.ABC):
    @abc.abstractmethod
    def get_function_space(self, mesh: fenics.Mesh) -> fenics.FunctionSpace:
        pass