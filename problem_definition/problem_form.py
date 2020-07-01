import abc

class ProblemForm(abc.ABC):
    def __init__(self):
        pass

    def k_tilde(self):
        pass

    def mass_form_m(self):
        pass

    def stiffness_form_k(self):
        pass

    def damping_form_c(self):
        pass

    def get_weak_form(self):
        pass