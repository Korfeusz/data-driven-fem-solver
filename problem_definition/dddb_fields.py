from typing import Optional

import fenics
import numpy as np
from space_definition import Spaces
from .fields import Fields

class DDDbFields(Fields):
    def __init__(self):
        super().__init__()
        self._old_constitutive_relation_multiplicative_parameters: Optional[fenics.Function] = None
        self._new_constitutive_relation_multiplicative_parameters: Optional[fenics.Function] = None
        self.tensor_space = None
        self._imported_displacement_field: Optional[fenics.Function] = None

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

    @new_constitutive_relation_multiplicative_parameters.setter
    def new_constitutive_relation_multiplicative_parameters(self, array: np.ndarray):
        self._new_constitutive_relation_multiplicative_parameters.vector()[:] = array

    @property
    def imported_displacement_field(self) -> fenics.Function:
        if self._imported_displacement_field is None:
            try:
                self._imported_displacement_field = fenics.Function(self.vector_space, name='imported_displacement')
            except TypeError:
                print('V space was not initialized before access')
        return self._imported_displacement_field


    def initialize(self, spaces: Spaces) -> None:
        self.vector_space = spaces.vector_space
        self.tensor_space = spaces.tensor_space