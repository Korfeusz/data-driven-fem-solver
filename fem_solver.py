from abc import ABC, abstractmethod

class FEMSolver(ABC):
    def __init__(self, mesh):
        self.mesh = mesh

    @abstractmethod
    def define_boundary_conditions(self):
        pass

