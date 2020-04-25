from abc import ABC, abstractmethod

class FEMSolver(ABC):
    def __init__(self, mesh_loader):
        self.mesh_loader = mesh_loader

    @abstractmethod
    def define_boundary_conditions(self):
        pass

