import abc

class Fields(abc.ABC):
    @abc.abstractmethod
    def generate(self, *args, **kwargs):
        pass