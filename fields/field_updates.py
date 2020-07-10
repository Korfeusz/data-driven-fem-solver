import abc

from .fields import Fields


class FieldUpdates(abc.ABC):
    @abc.abstractmethod
    def run(self, fields: Fields) -> None:
        pass