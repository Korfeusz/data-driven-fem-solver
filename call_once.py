class CallOnce:
    def __init__(self, fun):
        self.fun = fun
        self._called = False
        self._value = None

    def __call__(self, *args, **kwargs):
        if not self._called:
            self._called = True
            self._value = self.fun(*args, **kwargs)
        return self._value