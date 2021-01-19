from __future__ import print_function
from fenics import *
from dolfin_adjoint import *

# Define discrete Functionspace
mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Functions
u = TrialFunction(V)
v = TestFunction(V)
s = Function(V)                                          # PDE solution
lmbd = Function(V)                                       # Adjoint PDE solution
f = Function(V)                                          # Control parameter
alpha = Constant(1e-6)                                   # Regularisation parameter
ud = Expression("sin(pi*x[0])*sin(pi*x[1])", degree=4)   # Desired temperature profile

# Define variational problem
a = inner(grad(u), grad(v))*dx
L = f*v*dx

bcs = DirichletBC(V, 0.0, "on_boundary")

# Solve Poisson problem
solve(a == L, s, bcs)

J = Functional(0.5*(s-ud)**2*dx + alpha*f**2*dx)
m = Control(f)

rf = ReducedFunctional(J, m)

print("Difference between desired and actual heat profile before optimisation: %2.f" % errornorm(ud, s))

m_opt = minimize(rf, method="L-BFGS-B", tol=1e-10)
f.assign(m_opt)
solve(a == L, s, bcs)

print("Difference between desired and actual heat profile after optimisation: %.2f" % errornorm(ud, s))

plot(m_opt, title="Optimised control")
plot(s, title="Optimised profile")
plot(ud, title="Desired profile", mesh=mesh)
interactive()
