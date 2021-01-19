from fenics import *
from dolfin_adjoint import *

# Define discrete Functionspace
mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Functions
u = TrialFunction(V)
v = TestFunction(V)
s = Function(V)                   # PDE solution
lmbd = Function(V)                # Adjoint PDE solution
f = Function(V)                   # Control
alpha = Constant(1e-6)            # Regularisation parameter
ud = Expression("sin(pi*x[0])*sin(pi*x[1])", degree=4)   # Desired temperature profile

# Define variational problem
a = inner(grad(u), grad(v))*dx
L = f*v*dx

bcs = DirichletBC(V, 0.0, "on_boundary")

# Solve Poisson problem
solve(a == L, s, bcs)

plot(s, title="Temperature")
interactive()
