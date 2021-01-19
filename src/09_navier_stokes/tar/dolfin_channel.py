from __future__ import print_function
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-17"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

# Modified by Marie E. Rognes (meg@simula.no)

# Last changed: 2016-06-22

from fenics import *

# Adjust log level
set_log_level(PROGRESS)

# Parameters for time-stepping
dt = 0.0005
T = 0.1

# Material parameters
rho = 1000.0   # kg / m^3   density of water
mu  = 0.001002 # kg / (m*s) viscosity of water

# Load mesh from file
mesh = Mesh("dolfin_channel.xml.gz")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define boundary for inflow
def inflow_domain(x):
    return x[0] < DOLFIN_EPS and x[1] > DOLFIN_EPS and x[1] < (1 - DOLFIN_EPS)

# Define boundary for outflow
def outflow_domain(x):
    return x[0] > 1 - DOLFIN_EPS and x[1] > DOLFIN_EPS and x[1] < 1 - DOLFIN_EPS

# Define boundary for no-slip boundary condition
def noslip_domain(x, on_boundary):
    return on_boundary and not inflow_domain(x) and not outflow_domain(x)

# Define boundary conditions
inflow_bc  = DirichletBC(Q, 1e3, inflow_domain)
outflow_bc = DirichletBC(Q, 0, outflow_domain)
noslip_bc  = DirichletBC(V, (0, 0), noslip_domain)

# Collect boundary conditions
velocity_bcs = [noslip_bc]
pressure_bcs = [inflow_bc, outflow_bc]

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Define functions
u0 = Function(V)
u1 = Function(V)
p0 = Function(Q)
p1 = Function(Q)
U  = 0.5*(u0 + u)

# Define constants
rho = Constant(rho)
mu  = Constant(mu)
k   = Constant(dt)
f   = Constant((0, 0))

# Define facet normal
n = FacetNormal(mesh)

# Define Cauchy stress tensor
def sigma(u, p):
    return 2.0*mu*sym(grad(u)) - p*Identity(len(u))

# Tentative velocity step
F1 = rho*dot((u - u0)/k, v)*dx + rho*dot(grad(u0)*u0, v)*dx \
    + inner(sigma(U, p0), sym(grad(v)))*dx - dot(f, v)*dx \
    + dot(p0*n, v)*ds - mu*inner(grad(U).T*n, v)*ds
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure correction
a2 = dot(k*grad(p),  grad(q))*dx
L2 = dot(k*grad(p0), grad(q))*dx - rho*div(u1)*q*dx

# Velocity correction
a3 = dot(v, rho*u)*dx
L3 = dot(v, rho*u1)*dx + dot(v, k*grad(p0 - p1))*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Create files for storing solution
velocity_file = File("results/velocity.pvd")
pressure_file = File("results/pressure.pvd")

# Time-stepping
t = dt
progress = Progress("Time-stepping")
while t < T + DOLFIN_EPS:

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in velocity_bcs]
    solve(A1, u1.vector(), b1, "gmres", "ilu")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in pressure_bcs]
    solve(A2, p1.vector(), b2, "gmres", "amg")
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in velocity_bcs]
    solve(A3, u1.vector(), b3, "cg", "default")
    end()

    # Save solution to file
    velocity_file << u1
    pressure_file << p1

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u1, title="Velocity", rescale=True)

    # Move to next interval
    t += dt
    u0.assign(u1)
    p0.assign(p1)
    progress.update(t / T)

# Compute average velocity in x-direction
velocity = assemble(u1[0]*dx) / assemble(1.0*dx(mesh))
print("\nAverage x-velocity: %.16g" % velocity)

# Hold plot
interactive()
