import fenics

from init_safe import InitSafe
from .normal_function_space_creator import NormalFunctionSpaceCreator
from .vector_function_space_creator import VectorFunctionSpaceCreator
from .tensor_function_space_creator import TensorFunctionSpaceCreator

class Spaces(InitSafe):
    def __init__(self, vector_space_creator: VectorFunctionSpaceCreator,
                 tensor_space_creator: TensorFunctionSpaceCreator,
                 function_space_creator: NormalFunctionSpaceCreator = None):
        super().__init__()
        self.vector_space_creator = vector_space_creator
        self.tensor_space_creator = tensor_space_creator
        self.function_space_creator = function_space_creator
        self.vector_space = None
        self.tensor_space = None
        self.function_space = None

    def initialize(self, mesh: fenics.Mesh) -> None:
        self.vector_space = self.vector_space_creator.get_function_space(mesh=mesh)
        self.tensor_space = self.tensor_space_creator.get_function_space(mesh=mesh)
        if self.function_space_creator is not None:
            self.function_space = self.function_space_creator.get_function_space(mesh=mesh)


