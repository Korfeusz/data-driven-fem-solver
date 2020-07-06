import numpy as np


class NPYFile:
    def __init__(self, file_name: str):
        self.file_name = file_name
        self.file = open(file=file_name, mode='ab')

    def write(self, array: np.ndarray):
        np.save(self.file, array)


    def load(self, iteration):
        pass

    def close(self):
        self.file.close()
