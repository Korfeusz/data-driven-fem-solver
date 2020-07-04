from .vector_function_space_creator import VectorFunctionSpaceCreator
from .tensor_function_space_creator import TensorFunctionSpaceCreator

class Spaces:
    def __init__(self, vector_space_creator: VectorFunctionSpaceCreator,
                 tensor_space_creator: TensorFunctionSpaceCreator):
        self.vector_space_creator = vector_space_creator
        self.tensor_space_creator = tensor_space_creator
        self.vector_space = None
        self.tensor_space = None

    def generate(self, mesh):
        self.vector_space = self.vector_space_creator.get_function_space(mesh=mesh)
        self.tensor_space = self.tensor_space_creator.get_function_space(mesh=mesh)
