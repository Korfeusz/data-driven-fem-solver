import abc

class FieldUpdates(abc.ABC):
    @abc.abstractmethod
    def update_fields(self):
        pass