import fenics
from .marker_names import MarkerNames

IDLE_MARKER = 9999

class BoundaryMarkers:
    def __init__(self, tolerance):
        self.value = []
        self.tolerance = tolerance

    def mark_boundaries(self, boundary_definition: fenics.SubDomain, boundary_definition_settings: dict,
                        boundary_names: MarkerNames, mesh: fenics.mesh):
        self.value = fenics.MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
        self.value.set_all(IDLE_MARKER)
        for boundary_name in boundary_names:
            boundary = boundary_definition(self.tolerance, boundary_definition_settings[boundary_name])
            boundary.mark(self.value, boundary_name.value)
