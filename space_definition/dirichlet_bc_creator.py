from typing import List, Union
import fenics
from .marker_names import MarkerNames

class DirichletBCCreator:
    def __init__(self, boundary_names_for_condition: Union[List[MarkerNames], MarkerNames], value: fenics.Constant):
        if issubclass(type(boundary_names_for_condition), MarkerNames):
            self.boundary_names = [boundary_names_for_condition]
        else:
            self.boundary_names = boundary_names_for_condition
        self.value = value

    def apply(self, vector_space: fenics.VectorFunctionSpace, boundary_markers: fenics.MeshFunction):
        bcs = []
        for name in self.boundary_names:
            bcs.append(fenics.DirichletBC(vector_space, self.value, boundary_markers, name.value))
        return bcs
