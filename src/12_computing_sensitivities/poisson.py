from __future__ import print_function
from fenics import *

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

print("J = %.2f" % assemble((s-ud)**2*dx))
