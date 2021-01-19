from fenics import *
from dolfin_adjoint import *

# Define discrete Functionspace
mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Functions
u = TrialFunction(V)
v = TestFunction(V)
s = Function(V)                             # PDE solution
lmbd = Function(V)                          # Adjoint PDE solution
m = interpolate(Constant(1), V)             # Control parameter
ud = Expression("sin(pi*x[0])", degree=4)   # Desired temperature profile

# Define variational problem
a = inner(grad(u), grad(v))*dx
L = m*v*dx

bcs = DirichletBC(V, 0.0, "on_boundary")

# Solve Poisson problem
solve(a == L, s, bcs)

J = Functional(0.5*(s-ud)**2*dx)
m = Control(m)
dJdm = compute_gradient(J, m, project=True)

plot(dJdm)
interactive()
