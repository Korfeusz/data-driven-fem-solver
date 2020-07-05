import fenics

from space_definition import Spaces
from .fields import Fields

class DDDbFields(Fields):
    def __init__(self):
        super().__init__()
        self._old_constitutive_relation_multiplicative_parameters = None
        self._new_constitutive_relation_multiplicative_parameters = None
        self.tensor_space = None

    @property
    def old_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        if self._old_constitutive_relation_multiplicative_parameters is None:
            try:
                self._old_constitutive_relation_multiplicative_parameters = fenics.Function(self.tensor_space,
                                                                                            name='old_db_params')
            except TypeError:
                print('Tensor space was not initialized before access')
        return self._old_constitutive_relation_multiplicative_parameters

    @property
    def new_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        if self._new_constitutive_relation_multiplicative_parameters is None:
            try:
                self._new_constitutive_relation_multiplicative_parameters = fenics.Function(self.tensor_space,
                                                                                            name='new_db_params')
            except TypeError:
                print('Tensor space was not initialized before access')
        return self._new_constitutive_relation_multiplicative_parameters

    def initialize(self, spaces: Spaces) -> None:
        self.vector_space = spaces.vector_space
        self.tensor_space = spaces.tensor_space