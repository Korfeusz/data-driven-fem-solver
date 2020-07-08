import abc


class InitSafe(abc.ABC):
    def __init__(self):
        self._initialized = False

    def __getattribute__(self, item):
        if object.__getattribute__(self, '_initialized') or item == 'initialize':
            object.__setattr__(self, '_initialized', True)
        else:
            raise RuntimeError('Class {} has not been initialized before trying to access parameter {}'.format(type(self).__name__, item))
        return object.__getattribute__(self, item)

    @abc.abstractmethod
    def initialize(self, *args, **kwargs):
        pass