import abc
import fenics

from init_safe import InitSafe
from space_definition import Spaces

class Fields(InitSafe):
    @abc.abstractmethod
    def __init__(self):
        super().__init__()
        self._w = None
        self._u = None
        self._u_old = None
        self._v_old = None
        self._a_old = None
        self._u_new = None
        self.vector_space = None

    @abc.abstractmethod
    def initialize(self, spaces: Spaces) -> None:
        pass

    @property
    def w(self) -> fenics.Function:
        if self._w is None:
            try:
                self._w =  fenics.TestFunction(self.vector_space)
            except TypeError:
                print('Vector space was not initialized before access')
        return self._w

    @property
    def u(self) -> fenics.Function:
        if self._u is None:
            try:
                self._u =  fenics.TrialFunction(self.vector_space)
            except TypeError:
                print('Vector space was not initialized before access')
        return self._u

    @property
    def u_old(self) -> fenics.Function:
        if self._u_old is None:
            try:
                self._u_old =  fenics.Function(self.vector_space, name='u_old')
            except TypeError:
                print('Vector space was not initialized before access')
        return self._u_old

    @property
    def v_old(self) -> fenics.Function:
        if self._v_old is None:
            try:
                self._v_old =  fenics.Function(self.vector_space, name='v_old')
            except TypeError:
                print('Vector space was not initialized before access')
        return self._v_old

    @property
    def a_old(self) -> fenics.Function:
        if self._a_old is None:
            try:
                self._a_old  = fenics.Function(self.vector_space, name='a_old')
            except TypeError:
                print('Vector space was not initialized before access')
        return self._a_old

    @property
    def u_new(self) -> fenics.Function:
        if self._u_new is None:
            try:
                self._u_new =  fenics.Function(self.vector_space, name='displacement')
            except TypeError:
                print('Vector space was not initialized before access')
        return self._u_new
