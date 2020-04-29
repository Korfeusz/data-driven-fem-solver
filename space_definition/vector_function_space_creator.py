import fenics

class VectorFunctionSpaceCreator:
    def get_function_space(self, mesh: fenics.Mesh) -> fenics.FunctionSpace:
        return fenics.VectorFunctionSpace(mesh, 'P', 1)