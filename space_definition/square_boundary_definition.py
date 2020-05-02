import fenics
from .marker_index import MarkerIndex

class SquareSideMarkerIndex(MarkerIndex):
    X0 = ()
    X1 = ()
    Y0 = ()
    Y1 = ()

square_side_locations = {SquareSideMarkerIndex.X0: (0, 0),
                         SquareSideMarkerIndex.X1: (0, 1),
                         SquareSideMarkerIndex.Y0: (1, 0),
                         SquareSideMarkerIndex.Y1: (1, 1)}

class SquareBoundaryDefinition(fenics.SubDomain):
    def __init__(self, tolerance: float, side: tuple):
        super().__init__()
        self.tolerance = tolerance
        self.x_or_y = side[0]
        self.start_or_end = side[1]

    def inside(self, x, on_boundary):
        return on_boundary and fenics.near(x[self.x_or_y], self.start_or_end, self.tolerance)


def get_all_square_boundaries_and_marker_values():
    return {x.name: {'boundary': SquareBoundaryDefinition(tolerance=1e-14, side=x.value), 'marker_value': i} for i, x in enumerate(SquareSide)}
