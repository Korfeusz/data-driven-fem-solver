from typing import Optional

import fenics
import numpy as np

from lazy import Lazy
from space_definition import Spaces
from .fields import Fields
from .data_driven_parameters_space_names import DataDrivenParametersSpaceName

class DDDbFields(Fields):
    def __init__(self, parameters_space_name: DataDrivenParametersSpaceName):
        super().__init__()
        self._old_constitutive_relation_multiplicative_parameters: Optional[fenics.Function] = None
        self._new_constitutive_relation_multiplicative_parameters: Optional[fenics.Function] = None
        self.tensor_space = None
        self.function_space = None
        self._imported_displacement_field: Optional[fenics.Function] = None
        self._parameters_space_name = parameters_space_name


    def get_parameter_space(self) -> Optional[fenics.FunctionSpace]:
        space = None
        if self._parameters_space_name == DataDrivenParametersSpaceName.FUNCTION_SPACE:
            space = self.function_space
        elif self._parameters_space_name == DataDrivenParametersSpaceName.VECTOR_SPACE:
            space = self.vector_space
        elif self._parameters_space_name == DataDrivenParametersSpaceName.TENSOR_SPACE:
            space = self.tensor_space
        return space

    @property
    @Lazy
    def old_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        return fenics.Function(self.get_parameter_space(), name='old_db_params')


    @property
    def new_constitutive_relation_multiplicative_parameters(self) -> fenics.Function:
        if self._new_constitutive_relation_multiplicative_parameters is None:
            self._new_constitutive_relation_multiplicative_parameters = fenics.Function(self.get_parameter_space(),
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
        self.function_space = spaces.function_space