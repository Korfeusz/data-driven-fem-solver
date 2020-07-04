from .fem_solver import FemSolver
import fenics


class PreAssembledSolver(FemSolver):
    def __init__(self, problem, fields, boundary_conditions):
        super().__init__(problem, fields, boundary_conditions)
        self.a_form = problem.get_weak_form_lhs(fields=fields)
        self.L_form = problem.get_weak_form_rhs(fields=fields)
        self.boundary_conditions = boundary_conditions
        self.K, _ = fenics.assemble_system(self.a_form, self.L_form, self.boundary_conditions)
        self.solver = fenics.LUSolver(self.K, "mumps")
        self.solver.parameters["symmetric"] = True


    def run(self, fields):
        res = fenics.assemble(self.L_form)
        self.boundary_conditions[0].apply(res)
        self.solver.solve(self.K, fields.u_new.vector(), res)

