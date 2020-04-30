import fenics
from enum import Enum

class SquareSide(Enum):
    X0 = (0, 0)
    X1 = (0, 1)
    Y0 = (1, 0)
    Y1 = (1, 1)



class SquareBoundaryDefiner(fenics.SubDomain):
    def __init__(self, tolerance: float, side: SquareSide):
        super().__init__()
        self.tolerance = tolerance
        self.x_or_y = side.value[0]
        self.start_or_end = side.value[1]

    def inside(self, x, on_boundary):
        return on_boundary and fenics.near(x[self.x_or_y], self.start_or_end, self.tolerance)
