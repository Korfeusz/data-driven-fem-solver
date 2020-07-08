import abc


class InitSafe(abc.ABC):
    def __init__(self):
        self.initialized = False

    def __getattribute__(self, item):
        if object.__getattribute__(self, 'initialized') or item == 'initialize':
            object.__setattr__(self, 'initialized', True)
        else:
            raise RuntimeError('Class {} has not been initialized before trying to access parameter {}'.format(type(self).__name__, item))
        return object.__getattribute__(self, item)

    @abc.abstractmethod
    def initialize(self):
        pass