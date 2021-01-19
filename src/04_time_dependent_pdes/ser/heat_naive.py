from fenics import *

# Mesh and function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "Lagrange", 1)

# Time variables
dt = Constant(0.3)
t = float(dt)
T = 1.8

# Analytical solution
g_expr = "1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t"
g = Expression(g_expr, alpha=3.0, beta=1.2, t=0,
               degree=2)

# Define boundary condition
bc = DirichletBC(V, g, "on_boundary")

# Previous and current solution
u0 = interpolate(g, V)
u1 = Function(V)

# Define some functions
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.2 - 2.0 - 2*3.0)

# Variational problem at each time
a = u*v*dx + dt*inner(grad(u), grad(v))*dx
L = u0*v*dx + dt*f*v*dx

while t <= T:
    # Solve
    g.t = t  # Set g.t before it's evaluation inside solve
    solve(a == L, u1, bc)

    # Update for next timestep
    u0.assign(u1)
    t += float(dt)

plot(u1, title="Approximated final solution")
interactive()
