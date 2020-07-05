import fenics

from .constitutive_relation import ConstitutiveRelation

class DDDbConstitutiveRelation(ConstitutiveRelation):

    def get_new_value(self, r: fenics.Function) -> fenics.Expression:
        pass


    def get_old_value(self, r: fenics.Function) -> fenics.Expression:
        pass