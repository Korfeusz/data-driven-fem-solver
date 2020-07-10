import abc
from fields import Fields

class ProblemForm(abc.ABC):
    @abc.abstractmethod
    def get_weak_form_lhs(self, fields: Fields):
        pass

    @abc.abstractmethod
    def get_weak_form_rhs(self, fields: Fields):
        pass