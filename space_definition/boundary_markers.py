import fenics
from .marker_names import MarkerNames
from typing import Type

IDLE_MARKER = 9999

class BoundaryMarkers:
    def __init__(self, tolerance, boundary_definition_constructor: Type[fenics.SubDomain], boundary_definition_settings: dict,
                 boundary_names: Type[MarkerNames]):
        self.value = []
        self.tolerance = tolerance
        self.boundary_definition_constructor = boundary_definition_constructor
        self.boundary_definition_settings = boundary_definition_settings
        self.boundary_names = boundary_names

    def mark_boundaries(self, mesh: fenics.mesh):
        self.value = fenics.MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
        self.value.set_all(IDLE_MARKER)
        for boundary_name in self.boundary_names:
            boundary = self.boundary_definition_constructor(self.tolerance, self.boundary_definition_settings[boundary_name])
            boundary.mark(self.value, boundary_name.value)
