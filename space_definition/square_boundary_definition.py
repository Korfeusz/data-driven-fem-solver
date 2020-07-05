import fenics
from .marker_names import MarkerNames

class SquareSideMarkerNames(MarkerNames):
    X0 = ()
    X1 = ()
    Y0 = ()
    Y1 = ()

square_side_settings = {SquareSideMarkerNames.X0: (0, 0),
                        SquareSideMarkerNames.X1: (0, 1),
                        SquareSideMarkerNames.Y0: (1, 0),
                        SquareSideMarkerNames.Y1: (1, 1)}

class SquareBoundaryDefinition(fenics.SubDomain):
    def __init__(self, tolerance: float, side: tuple):
        super().__init__()
        self.tolerance = tolerance
        self.x_or_y = side[0]
        self.start_or_end = side[1]

    def inside(self, x, on_boundary):
        return on_boundary and fenics.near(x[self.x_or_y], self.start_or_end, self.tolerance)
