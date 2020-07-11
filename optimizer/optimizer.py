import abc

class Optimizer(abc.ABC):
    @abc.abstractmethod
    def function_to_minimize(self, parameters):
        pass

    @abc.abstractmethod
    def run(self):
        pass
