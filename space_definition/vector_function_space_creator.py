import fenics
from .function_space_creator import FunctionSpaceCreator

class VectorFunctionSpaceCreator(FunctionSpaceCreator):
    def __init__(self, element_family: str, degree: int):
        self.element_family = element_family
        self.degree = degree

    def get_function_space(self, mesh: fenics.Mesh) -> fenics.FunctionSpace:
        return fenics.VectorFunctionSpace(mesh, self.element_family, self.degree)