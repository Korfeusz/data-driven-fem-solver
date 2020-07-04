import abc



class FieldUpdates(abc.ABC):
    @abc.abstractmethod
    def run(self, *args, **kwargs):
        pass