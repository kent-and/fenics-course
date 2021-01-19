from fenics import *
from dolfin_adjoint import *

# Define mesh and finite element space
mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define basis functions and parameters
u = TrialFunction(V)
v = TestFunction(V)
m = interpolate(Constant(1.0), V)
nu = Constant(1.0)

# Define variational problem
a = nu*inner(grad(u), grad(v))*dx
L = m*v*dx
bc = DirichletBC(V, 0.0, "on_boundary")

# Solve variational problem
u = Function(V)
solve(a == L, u, bc)
plot(u, title="u")

# Define 'observations'
u_d = Expression("sin(pi*x[0])*sin(pi*x[1])", degree=4)

# Define 'misfit' functional
j = 0.5*(u - u_d)*(u - u_d)*dx
J = Functional(j)

# Define the control variable
m = Control(m)

# Compute gradient of J with respect to m (dJ/dm):
dJdm = compute_gradient(J, m, project=True, forget=False)
plot(dJdm, title="dJdm", interactive=True)

# Run optimization
R = ReducedFunctional(J, m)
m_opt = minimize(R)
plot(m_opt, title="Optimal control", interactive=True)

# Compute Hessian:
H = hessian(J, m)
direction = interpolate(Constant(1.0), V)
plot(H(direction), title="H( 1 )")
interactive()
