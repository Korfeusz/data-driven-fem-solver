import fenics


class HDF5File:
    def __init__(self, mesh: fenics.Mesh, mode: str, file_name: str, function: fenics.Function, function_name: str):
        self.file = fenics.HDF5File(mesh.mpi_comm(), file_name, mode)
        self.function = function
        self.function_name = function_name



    def save(self,  iteration: int):
        self.file.write(self.function, self.function_name, iteration)


    def load(self, iteration):
        self.file.read(self.function, '{}/vector_{}'.format(self.function_name, iteration))