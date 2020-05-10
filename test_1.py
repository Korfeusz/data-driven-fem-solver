from space_definition import UnitSquareMeshCreator, MeshCreator, FunctionSpaceCreator, \
    VectorFunctionSpaceCreator, SquareBoundaryDefinition, BoundaryMarkers, square_side_settings,\
    SquareSideMarkerNames, DirichletBCCreator, TensorFunctionSpaceCreator
import fenics
import numpy as np

mesh_creator = UnitSquareMeshCreator(horizontal_cell_no=5, vertical_cell_no=5)
vector_space_creator = VectorFunctionSpaceCreator(element_family='CG', degree=1)
tensor_space_creator = TensorFunctionSpaceCreator(element_family='DG', degree=0)
tolerance = 1e-14

boundary_markers = BoundaryMarkers(tolerance=tolerance, boundary_definition_constructor=SquareBoundaryDefinition,
                                   boundary_definition_settings=square_side_settings,
                                   boundary_names=SquareSideMarkerNames)
bc_creator = DirichletBCCreator(boundary_names_for_condition=SquareSideMarkerNames.X0,
                                value=fenics.Constant((0, 0)))


mesh = mesh_creator.get_mesh()
vector_space = vector_space_creator.get_function_space(mesh=mesh)
tensor_space = tensor_space_creator.get_function_space(mesh=mesh)
boundary_markers.mark_boundaries(mesh=mesh)
bc = bc_creator.apply(vector_space=vector_space, boundary_markers=boundary_markers.value)
ds = fenics.Measure('ds', domain=mesh, subdomain_data=boundary_markers.value)

# Elastic parameters
E  = 1000.0
nu = 0.3
mu    = fenics.Constant(E / (2.0*(1.0 + nu)))
lmbda = fenics.Constant(E*nu / ((1.0 + nu)*(1.0 - 2.0*nu)))

# Mass density
rho = fenics.Constant(1.0)

# Generalized-alpha method parameters
alpha_m = fenics.Constant(0.2)
alpha_f = fenics.Constant(0.4)
gamma   = fenics.Constant(0.5+alpha_f-alpha_m)
beta    = fenics.Constant((gamma+0.5)**2/4.)


# Time-stepping parameters
T       = 4.0
Nsteps  = 50
dt = fenics.Constant(T/Nsteps)
# Time-stepping
time = np.linspace(0, T, Nsteps+1)


p0 = 1.
cutoff_Tc = T/5
# Define the loading as an expression depending on t
p = fenics.Expression(("t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=1)

# Test and trial functions
du = fenics.TrialFunction(vector_space)
u_ = fenics.TestFunction(vector_space)
# Current (unknown) displacement
u = fenics.Function(vector_space, name="Displacement")
# Fields from previous time step (displacement, velocity, acceleration)
u_old = fenics.Function(vector_space)
v_old = fenics.Function(vector_space)
a_old = fenics.Function(vector_space)


# Stress tensor
def sigma(r):
    return 2.0*mu*fenics.sym(fenics.grad(r)) + lmbda*fenics.tr(fenics.sym(fenics.grad(r)))*fenics.Identity(len(r))

def np_sym(r):
    return 0.5 * (r + np.transpose(r))

def sigma_np(r, mu, lmbda):
    return 2.0 * mu * np_sym(np.gradient(r)) + lmbda * np.trace(np_sym(np.gradient(r))) * np.identity(len(r))


s = fenics.Function(tensor_space)
vv = fenics.FunctionSpace(mesh, 'CG', 1)
class Sigma(fenics.UserExpression):
    def set_u(self, r):
        self.r = r


    def eval(self, value, x):
        print('x ', x)
        print('val ', value)
        # print('s ', fenics.project(self.r, V=vector_space))
        # print('r', self.r[x[0]][x[1]])
        # print('type: ', type(sigma(self.r)))
        value[0], value[1], value[2], value[3] = 0.1, 0.1, 0.1, 0.1

    def value_shape(self):
        return (2, 2)

sig = Sigma()
# Mass form
def m(u, u_):
    return rho*fenics.inner(u, u_)*fenics.dx

# Elastic stiffness form
# def k(u, u_):
#     return fenics.inner(sigma(u), fenics.sym(fenics.grad(u_)))*fenics.dx

def k(u, u_):
    print(type(u))
    sig.set_u(u)
    return fenics.inner(sig, fenics.sym(fenics.grad(u_)))*fenics.dx


# Work of external forces
def Wext(u_):
    return fenics.dot(u_, p)*ds(SquareSideMarkerNames.X1.value)

# Update formula for acceleration
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return (u-u_old-dt_*v_old)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a_old

# Update formula for velocity
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u, u_old, v_old, a_old):
    """Update fields at the end of each time step."""

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # use update functions using vector arguments
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # Update (u_old <- u)
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()

def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new


#
# # Define solver for reusing factorization
# K, res = fenics.assemble_system(a_form, L_form, bc)
# solver = fenics.LUSolver(K, "mumps")
# solver.parameters["symmetric"] = True



xdmf_file = fenics.XDMFFile("elastodynamics-results.xdmf")
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = False


for (i, dt) in enumerate(np.diff(time)):

    t = time[i+1]
    print("Time: ", t)

    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    p.t = t-float(alpha_f*dt)

    # Residual
    a_new = update_a(du, u_old, v_old, a_old, ufl=False)
    v_new = update_v(a_new, u_old, v_old, a_old, ufl=False)
    res = m(avg(a_old, a_new, alpha_m), u_) + k(avg(u_old, du, alpha_f), u_) - Wext(u_)
    a_form = fenics.lhs(res)
    L_form = fenics.rhs(res)
    # Define solver for reusing factorization
    K, res = fenics.assemble_system(a_form, L_form, bc)
    solver = fenics.LUSolver(K, "mumps")
    solver.parameters["symmetric"] = True

    # Solve for new displacement
    # res = fenics.assemble(L_form)
    bc[0].apply(res)
    solver.solve(K, u.vector(), res)


    # Update old fields with new quantities
    update_fields(u, u_old, v_old, a_old)

    # Save solution to XDMF format
    xdmf_file.write(u, t)


    p.t = t