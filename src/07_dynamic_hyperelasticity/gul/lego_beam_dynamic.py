from __future__ import print_function
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-18"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

# Modified by Martin Alnaes
# Last changed: 2016-04-10

from fenics import *

# Adjust log level
set_log_level(PROGRESS)

# Turn on C++ optimization of generated code
parameters["form_compiler"]["cpp_optimize"] = True

# Select the newer code generation backend 'uflacs'
parameters["form_compiler"]["representation"] = "uflacs"

# Parameters for time-stepping
dt = 0.002
T = 0.05

# Material parameters
rho   = 1450.0       # kg / m^3 density of PVD plastic
mu    = 0.0023 * 1e9 # N  / m^2 Lame parameter mu for PVC plastic
lmbda = 0.0105 * 1e9 # N  / m^2 Lame parameter lambda for PVC plastic

# Load mesh from file
mesh = Mesh("lego_beam.xml")
plot(mesh)

# Define function space (P1-P1)
P1 = VectorElement("Lagrange", tetrahedron , 1)
VV = FunctionSpace(mesh, MixedElement([P1, P1]))

# Define left boundary
class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 0.1*0.001 + DOLFIN_EPS

# Define right boundary
class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 10.0*0.008 - 0.1*0.001 - DOLFIN_EPS

# Mark boundaries
left_boundary = LeftBoundary()
right_boundary = RightBoundary()
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left_boundary.mark(boundaries, 1)
right_boundary.mark(boundaries, 2)

# Redefine boundary measure
ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# Define boundary conditions
clamp_bc = DirichletBC(VV.sub(0), (0, 0, 0), boundaries, 1)
bcs = [clamp_bc]

# Define forces
B = Constant((0, 0, -9.81*rho))
G = Constant((0, 0, -5e3)) # using G in place of T since T is already used

# Define functions and test functions
up0 = Function(VV)
up1 = Function(VV)
u0, p0 = split(up0)
u1, p1 = split(up1)
v, q = TestFunctions(VV)

# Define midpoint values
u = 0.5*(u0 + u1)
p = 0.5*(p0 + p1)

# Define constants
rho   = Constant(rho)
mu    = Constant(mu)
lmbda = Constant(lmbda)
k     = Constant(dt)

# Define strain measures
I = Identity(len(u)) # the identity matrix
F = I + grad(u)      # the deformation gradient
C = F.T*F            # the right Cauchy-Green deformation tensor
E = 0.5*(C - I)      # the Green-Lagrange strain tensor

# Define strain energy density
E = variable(E)
W = lmbda/2*(tr(E)**2) + mu*tr(E*E)

# Define Piola-Kirchoff stress tensors
S = diff(W, E) # the second Piola-Kirchoff stress tensor
P = F*S        # the first Piola-Kirchoff stress tensor

# Define nonlinear problem for one time-step
R = rho*dot(p1 - p0, v)*dx + k*inner(P, grad(v))*dx \
    - k*dot(B, v)*dx - k*dot(G, v)*ds(2) \
    + dot(u1 - u0, q)*dx - k*dot(p, q)*dx

# Compute the derivative outside of the time loop
J = derivative(R, up1)

# Create files for storing solution
displacement_file = File("results/displacement.pvd")
velocity_file = File("results/velocity.pvd")

# Time-stepping
t = dt
progress = Progress("Time-stepping")
while t < T + DOLFIN_EPS:

    # Solve nonlinear problem for time step
    solve(R == 0, up1, bcs, J=J)

    # Save solution to file
    _u1, _p1 = up1.split() # this should be fixed in DOLFIN
    displacement_file << _u1
    velocity_file << _p1

    # Move to next interval
    t += dt
    up0.assign(up1)
    progress.update(t / T)

# Plot solution
plot(u1, title="Displacement", mode="displacement")
plot(p1, title="Velocity")

# Compute average displacement in z-direction
displacement = assemble(u1[2]*dx) / assemble(1.0*dx(mesh))
print("\nAverage z-displacement: %.16g" % displacement)

# This example highlights some improvements that could be made in DOLFIN:
# 1. The use of split(w) vs w.split() is very confusing
# 2. Repeated plotting of a function created using split(w) creates new windows
# 3. Functions created using w.split() cannot be used with derivative()
# 4. Functions created using split(w) cannot be used for plotting
