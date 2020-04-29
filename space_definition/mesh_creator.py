import abc

class MeshCreator(abc.ABC):
    @abc.abstractmethod
    def get_mesh(self):
        pass