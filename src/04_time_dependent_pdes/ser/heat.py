from __future__ import print_function
from fenics import *

# Create mesh and define function space
N = 8
mesh = UnitSquareMesh(N, N)
q = 1
V = FunctionSpace(mesh, 'Lagrange', q)

# Define initial condition expression g (will also be used as boundary
# condition analytical solution), and interpolate into initial function u0
alpha = 3.0
beta = 1.2
t = Constant(0.0)
g = Expression('x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
               alpha=alpha, beta=beta, t=t,
               degree=2)
u0 = interpolate(g, V)

# Define boundary condition
bc = DirichletBC(V, g, "on_boundary")

# Define timestep and end-time
dt = Constant(0.1)
T = 1.8

# Define some functions
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

# Define variational problem for each time-step
a = u*v*dx + dt*inner(grad(u), grad(v))*dx
L = (u0 + dt*f)*v*dx

# Assemble once before the time-stepping begins
A = assemble(a)

# Define function for unknown at this time step
u1 = Function(V)

# Run time-loop
t.assign(dt)

while float(t) <= (T + 100*DOLFIN_EPS):

    print("t = %g" %  t)
    # Assemble right-hand side vector
    b = assemble(L)

    # Apply boundary condition
    bc.apply(A, b)

    # Solve linear system of equations
    solve(A, u1.vector(), b)

    # Update time and previous function
    t.assign(float(t + dt))
    u0.assign(u1)

    # Plot solution
    plot(u1)

# Make sure that t == T for computing error norms
t.assign(T)

# Plot final solution and exact solution
plot(u1, title="Approximated solution")
plot(g, mesh=mesh, title="Exact solution")

# Compute the l2 error norm (three almost equivalent ways):
rise = 2
E = errornorm(g, u1, degree_rise=rise)
#E = sqrt(assemble((g-u1)**2*dx))
#E = sqrt( assemble( (g-u1)**2*dx(degree=2*(q + rise)) ) )

print("error = %.14f" % E)
interactive()
