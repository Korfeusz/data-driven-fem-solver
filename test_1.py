from space_definition import UnitSquareMeshCreator, MeshCreator, FunctionSpaceCreator, \
    VectorFunctionSpaceCreator, SquareBoundaryDefinition, BoundaryMarkers, square_side_settings, SquareSideMarkerNames
import fenics

mesh = UnitSquareMeshCreator(horizontal_cell_no=94, vertical_cell_no=94).get_mesh()
V = VectorFunctionSpaceCreator(element_family='CG', degree=1).get_function_space(mesh)
print(type(fenics.FunctionSpace)  == type(V))
print(isinstance(V, fenics.FunctionSpace))
tolerance = 1e-14

boundary_markers = BoundaryMarkers(tolerance=tolerance, boundary_definition_constructor=SquareBoundaryDefinition,
                                   boundary_definition_settings=square_side_settings,
                                   boundary_names=SquareSideMarkerNames, ).mark_boundaries( mesh=mesh)

