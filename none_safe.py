class NoneSafe:
    def __getattribute__(self, item):
        element = object.__getattribute__(self, item)
        if element is None:
            raise RuntimeError('Class {} has not been initialized before trying to access parameter {}'.format(type(self).__name__, item))
        return element