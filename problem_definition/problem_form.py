import abc

class ProblemForm(abc.ABC):
    @abc.abstractmethod
    def get_weak_form_lhs(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def get_weak_form_rhs(self, *args, **kwargs):
        pass