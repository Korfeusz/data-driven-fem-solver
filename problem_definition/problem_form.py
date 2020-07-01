import fenics



class ProblemForm:
    def __init__(self, rho, eta_m, eta_k, constitutive_relation, gamma, beta, alpha_f, alpha_m):
        self.rho = rho
        self.eta_m = eta_m
        self.eta_k = eta_k
        self.constitutive_relation = constitutive_relation
        self.gamma = gamma
        self.beta = beta
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m

    def mass_form_m(self, u, w):
        return self.rho*fenics.inner(u, w)*fenics.dx

    def damping_form_c(self, u, w):
        return self.eta_m*self.mass_form_m(u, w) + self.eta_k*self.stiffness_form_k(u, w)

    def stiffness_form_k(self, u, w):
        return fenics.inner(self.constitutive_relation(u), fenics.sym(fenics.grad(w))) * fenics.dx

    def get_weak_form_lhs(self, u, w):
        k = self.stiffness_form_k
        m = self.mass_form_m
        c = self.damping_form_c
        # k(u, w) +