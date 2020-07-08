from typing import Optional

import fenics
import numpy as np

from lazy import Lazy
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
    @Lazy
    def old_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        return fenics.Function(self.tensor_space, name='old_db_params')

    @property
    def new_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        if self._new_constitutive_relation_multiplicative_parameters is None:
            self._new_constitutive_relation_multiplicative_parameters = fenics.Function(self.tensor_space,
                                                                                            name='new_db_params')
        return self._new_constitutive_relation_multiplicative_parameters

    @new_constitutive_relation_multiplicative_parameters.setter
    def new_constitutive_relation_multiplicative_parameters(self, array: np.ndarray):
        self._new_constitutive_relation_multiplicative_parameters.vector()[:] = array

    @property
    @Lazy
    def imported_displacement_field(self) -> fenics.Function:
        return fenics.Function(self.vector_space, name='imported_displacement')


    def initialize(self, spaces: Spaces) -> None:
        self.vector_space = spaces.vector_space
        self.tensor_space = spaces.tensor_space