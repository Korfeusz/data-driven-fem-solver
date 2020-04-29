import fenics
from .mesh_creator import MeshCreator


class UnitSquareMeshCreator(MeshCreator):
    def __init__(self, horizontal_cell_no: int, vertical_cell_no: int):
        self.horizontal_cell_no = horizontal_cell_no
        self.vertical_cell_no = vertical_cell_no

    def get_mesh(self) -> fenics.Mesh:
        return fenics.UnitSquareMesh(nx=self.horizontal_cell_no, ny=self.vertical_cell_no)