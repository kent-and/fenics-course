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

# Material parameters
rho   = 1450.0       # kg / m^3 density of PVC plastic
mu    = 0.0023 * 1e9 # N  / m^2 Lame parameter mu for PVC plastic
lmbda = 0.0105 * 1e9 # N  / m^2 Lame parameter lambda for PVC plastic

# Create mesh
n = 8
mesh = BoxMesh(Point(0, 0, 0), Point(5, 1, 1), n, n, 5*n)

# Define function space (P1)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Define left boundary
class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

# Define right boundary
class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 5.0)

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
clamp_bc = DirichletBC(V, (0, 0, 0), boundaries, 1)
bcs = [clamp_bc]

# Define forces
B = Constant((0, 0, -9.81*rho))
T = Constant((0, 0, -5e3))

# Define solution u and test function v
u = Function(V)
v = TestFunction(V)

# Define constants
rho   = Constant(rho)
mu    = Constant(mu)
lmbda = Constant(lmbda)

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

# Define nonlinear problem
R = inner(P, grad(v))*dx - dot(B, v)*dx - dot(T, v)*ds(2)

# Alternative formulation as first variation of total energy
#Q = W*dx - dot(B, u)*dx - dot(T, u)*ds(2)
#R = derivative(Q, u, v)

# Note: We would normally use 'F' to denote the nonlinear residual
# form but we use 'R' here to avoid confusion with the deformation
# gradient.

# Compute solution
solve(R == 0, u, bcs)

# Compute average displacement in z-direction
displacement = assemble(u[2]*dx) / assemble(1.0*dx(mesh))
print("\nAverage z-displacement: %.16g" % displacement)

# Save solution to file
file = File("displacement.pvd")
file << u

# Plot solution
plot(u, title="Displacement", mode="displacement")
interactive()
