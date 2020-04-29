import abc
import fenics

class MeshCreator(abc.ABC):
    @abc.abstractmethod
    def get_mesh(self) -> fenics.Mesh:
        pass