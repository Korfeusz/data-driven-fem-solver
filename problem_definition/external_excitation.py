import abc

class ExternalExcitation(abc.ABC):
    @abc.abstractmethod
    def value(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def update(self, *args, **kwargs):
        pass