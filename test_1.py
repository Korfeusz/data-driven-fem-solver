from space_definition import UnitSquareMeshCreator, MeshCreator, FunctionSpaceCreator, \
    VectorFunctionSpaceCreator, SquareBoundaryDefinition, SquareSide
import fenics

mesh = UnitSquareMeshCreator(horizontal_cell_no=94, vertical_cell_no=94).get_mesh()
V = VectorFunctionSpaceCreator(element_family='CG', degree=1).get_function_space(mesh)
print(type(fenics.FunctionSpace)  == type(V))
print(isinstance(V, fenics.FunctionSpace))

