import fenics
from constitutive_relations import ConstitutiveRelation


class ProblemForm:
    def __init__(self, rho, eta_m, eta_k, constitutive_relation: ConstitutiveRelation, gamma, beta, alpha_f, alpha_m, delta_t, f_ext):
        self.rho = rho
        self.eta_m = eta_m
        self.eta_k = eta_k
        self.constitutive_relation = constitutive_relation
        self.gamma = gamma
        self.beta = beta
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m
        self.delta_t = delta_t
        self.f_ext = f_ext

    @property
    def c_1(self):
        return self.gamma * (1 - self.alpha_f) / (self.beta * self.delta_t)

    @property
    def c_2(self):
        return 1 - self.gamma/self.beta * (1 - self.alpha_f)

    @property
    def c_3(self):
        return self.delta_t * (1 - self.alpha_f) * (1 - (self.gamma / (2 * self.beta)))

    @property
    def m_1(self):
        return  (1 - self.alpha_m) / (self.delta_t*self.delta_t * self.beta)

    @property
    def m_2(self):
        return  (1 - self.alpha_m) / (self.delta_t * self.beta)

    @property
    def m_3(self):
        return  1 - (1 - self.alpha_m) / (2 * self.beta)


    def m(self, u, w):
        return self.rho*fenics.inner(u, w)*fenics.dx

    def c(self, u, w, constitutive_relation_function):
        return self.eta_m * self.m(u, w) + self.eta_k * self.k(u, w, constitutive_relation_function)

    def k(self, u, w, constitutive_relation_function):
        return fenics.inner(constitutive_relation_function(u), fenics.sym(fenics.grad(w))) * fenics.dx

    def get_weak_form_lhs(self, u, w):
        return (1 - self.alpha_f) * self.k(u, w, self.constitutive_relation.get_new_value) \
               + self.c_1 * self.c(u, w, self.constitutive_relation.get_new_value) + self.m_1 * self.m(u, w)

    def get_weak_form_rhs(self, u, v, a, w):
        old_relation = self.constitutive_relation.get_old_value
        return -self.alpha_f * self.k(u, w, old_relation) + self.f_ext(w) \
               + self.c_1 * self.c(u, w, old_relation) - self.c_2 * self.c(v, w, old_relation)\
               - self.c_3 * self.c(a, w, old_relation)\
               + self.m_1 * self.m(u, w) + self.m_2 * self.m(v, w) - self.m_3 * self.m(a, w)