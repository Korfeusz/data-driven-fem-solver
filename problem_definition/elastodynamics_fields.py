import fenics
from .fields import Fields

class ElastodynamicsFields(Fields):
    def __init__(self):
        self.w = None
        self.u = None
        self.u_old = None
        self.v_old = None
        self.a_old = None
        self.u_new = None

    def generate(self, vector_space):
        self.w = fenics.TestFunction(vector_space)
        self.u = fenics.TrialFunction(vector_space)
        self.u_old = fenics.Function(vector_space, name='u_old')
        self.v_old = fenics.Function(vector_space, name='v_old')
        self.a_old = fenics.Function(vector_space, name='a_old')
        self.u_new = fenics.Function(vector_space, name='displacement')