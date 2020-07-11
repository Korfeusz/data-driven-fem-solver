import fenics
from fields import DDDbFields, DataDrivenParametersSpaceName
from .constitutive_relation import ConstitutiveRelation

class DDDbConstitutiveRelation(ConstitutiveRelation):
    def __init__(self, fields: DDDbFields):
        self.fields = fields

    def get_new_value(self, r: fenics.Function) -> fenics.Expression:
        return fenics.elem_mult(self.fields.new_constitutive_relation_multiplicative_parameters,
                                fenics.sym(fenics.grad(r)))


    def get_old_value(self, r: fenics.Function) -> fenics.Expression:
        return fenics.elem_mult(self.fields.old_constitutive_relation_multiplicative_parameters,
                                fenics.sym(fenics.grad(r)))
