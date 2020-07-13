from file_handling import NPYFile
from lazy import Lazy
import numpy as np

class DDDbLoader:
    def __init__(self, database_file_name: str, strains_file_prefix: str, params_file_prefix: str,
                 iteration_number: int, path_to_files=None):
        self.database_file_name = database_file_name
        self.strains_file_prefix = strains_file_prefix
        self.params_file_prefix = params_file_prefix
        self.path_to_files = path_to_files
        self.file = NPYFile(database_file_name)
        self.iteration_number = iteration_number


    @property
    @Lazy
    def parameters_array(self) -> np.ndarray:
        first_array =  self.file.load(self.params_file_prefix, iteration=0)
        if len(first_array.shape) == 1:
            return self.load_function(self.params_file_prefix, first_array)
        elif len(first_array.shape) == 3:
            return  self.load_tensor(self.params_file_prefix, first_array)

    @property
    @Lazy
    def strains_array(self) -> np.ndarray:
        first_array =  self.file.load(self.strains_file_prefix, iteration=0)
        return self.load_tensor(self.strains_file_prefix, first_array)


    def load_tensor(self, prefix: str, array: np.ndarray) -> np.ndarray:
        for iteration in range(1, self.iteration_number):
            array = np.vstack((array, self.file.load(prefix, iteration)))
        return array

    def load_function(self, prefix: str, array: np.ndarray) -> np.ndarray:
        for iteration in range(1, self.iteration_number):
            array = np.append(array, self.file.load(prefix, iteration))
        return array
