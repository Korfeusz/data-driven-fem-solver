import fenics
import numpy as np

mesh = fenics.UnitSquareMesh(nx=5, ny=5)
# print(mesh.cells())
V = fenics.VectorFunctionSpace(mesh, 'CG', 1)
# tensor_space = fenics.TensorFunctionSpace(mesh, 'DG', 0)

def left(x, on_boundary):
    return fenics.near(x[0], 0.) and on_boundary

def right(x, on_boundary):
    return fenics.near(x[0], 1.) and on_boundary

boundary_subdomains = fenics.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
force_boundary = fenics.AutoSubDomain(right)
force_boundary.mark(boundary_subdomains, 2)

ds = fenics.ds(subdomain_data=boundary_subdomains)

zero = fenics.Constant((0.0, 0.0))
bc = fenics.DirichletBC(V, zero, left)

fs = fenics.FunctionSpace(mesh, 'Lagrange', 1)
young_modulus = fenics.Function(fs, name='young_modulus')

v2d = fenics.vertex_to_dof_map(fs)

young_modulus.vector()[v2d] = np.linspace(0, 10, 36)

print(young_modulus.vector()[:])
print(young_modulus((0.35, 0.11)))

# element = fs.element()
# dofmap = fs.dofmap()
# for cell in fenics.cells(mesh):
#     print(element.tabulate_dof_coordinates(cell))
#     print(dofmap.cell_dofs(cell.index()))

def sigma(u, poisson_coefficient, young_modulus):
    return 2.0*young_modulus/(2*(1 + poisson_coefficient))*fenics.sym(fenics.grad(u))\
           + young_modulus * poisson_coefficient /((1+poisson_coefficient)*(1-2*poisson_coefficient))\
           *fenics.tr(fenics.sym(fenics.grad(u)))*fenics.Identity(len(u))

u = fenics.TrialFunction(V)
v = fenics.TestFunction(V)

poisson_coefficient = 0.3
a = fenics.inner(sigma(u, poisson_coefficient, young_modulus), fenics.sym(fenics.grad(v))) * fenics.dx
f = fenics.Constant((0, 0))
traction = fenics.Constant((1, 0))
L = fenics.dot(f, v) * fenics.dx + fenics.dot(traction, v) * (ds(2))
# F = a - L

u = fenics.Function(V)
fenics.solve(a == L, u, bc)
# fenics.solve(F==0, u, bc)
file = fenics.File("test/output.pvd")
file << u
file << young_modulus

import replace_uint