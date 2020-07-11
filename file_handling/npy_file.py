import numpy as np


class NPYFile:
    def __init__(self, file_name: str):
        self.file_name = file_name
        self.file_name = file_name

    def write(self, prefix: str, iteration: int, array: np.ndarray):
        np.save('{}_{}_{}.npy'.format(prefix, self.file_name, iteration), array)


    def load(self, prefix: str, iteration: int):
        return np.load('{}_{}_{}.npy'.format(prefix, self.file_name, iteration))

