from enum import Enum

class MarkerIndex(Enum):
    def __repr__(self):
        return '{}'.format(self.name)

    def __str__(self):
        return '{}'.format(self.name)

    def __new__(cls):
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj