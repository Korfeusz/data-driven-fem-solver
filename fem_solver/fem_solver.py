import abc


class FemSolver(abc.ABC):
    @abc.abstractmethod
    def prepare(self):
        pass

    @abc.abstractmethod
    def run(self):
        pass