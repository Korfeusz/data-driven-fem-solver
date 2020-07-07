import abc


class Builder(abc.ABC):
    @abc.abstractmethod
    def __init__(self, class_type):
        self.class_type = class_type

    @abc.abstractmethod
    def set(self, **kwargs):
        pass

    @abc.abstractmethod
    def build(self):
        pass