import abc
import fenics
from typing import Type

class ConstitutiveRelation(abc.ABC):
    @abc.abstractmethod
    def get_new_value(self, r):
        pass

    @abc.abstractmethod
    def get_old_value(self, r):
        pass