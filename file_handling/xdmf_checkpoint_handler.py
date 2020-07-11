import fenics

class XDMFCheckpointHandler:
    def __init__(self, file_name: str, clear_at_start: bool, field: fenics.Function, field_name: str):
        self.clear = clear_at_start
        self.file = fenics.XDMFFile(file_name)
        self.field = field
        self.field_name = field_name

    def write(self, iteration: int) -> None:
        self.file.write_checkpoint(self.field, self.field_name, iteration, fenics.XDMFFile.Encoding.HDF5, self.clear)
        self.clear = False

#
# f_out.write_checkpoint(project(f, V), "f", 0.0, XDMFFile.Encoding.HDF5, False)

# f_out.write_checkpoint(project(f, V), "f", t, XDMFFile.Encoding.HDF5, True)
