import numpy as np


class NPYFile:
    def __init__(self, file_name):
        self.file_name = file_name

    def write(self, array: np.ndarray):
        np.save(self.file_name, array)


    def load(self, iteration):
        pass
