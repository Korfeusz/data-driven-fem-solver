import abc
import fenics

from lazy import Lazy
from init_safe import InitSafe
from space_definition import Spaces

class Fields(InitSafe):
    @abc.abstractmethod
    def __init__(self):
        super().__init__()
        self.vector_space = None

    @abc.abstractmethod
    def initialize(self, spaces: Spaces) -> None:
        pass

    @property
    @Lazy
    def w(self) -> fenics.Function:
        return fenics.TestFunction(self.vector_space)

    @property
    @Lazy
    def u(self) -> fenics.Function:
        return fenics.TrialFunction(self.vector_space)

    @property
    @Lazy
    def u_old(self) -> fenics.Function:
        return fenics.Function(self.vector_space, name='u_old')

    @property
    @Lazy
    def v_old(self) -> fenics.Function:
        return fenics.Function(self.vector_space, name='v_old')

    @property
    @Lazy
    def a_old(self) -> fenics.Function:
        return fenics.Function(self.vector_space, name='a_old')

    @property
    @Lazy
    def u_new(self) -> fenics.Function:
        return fenics.Function(self.vector_space, name='displacement')
