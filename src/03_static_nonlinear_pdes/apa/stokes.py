from __future__ import print_function
from fenics import *

# Use -02 optimization
parameters["form_compiler"]["cpp_optimize"] = True

# Define mesh and geometry
mesh = Mesh("dolfin-2.xml.gz")
x = SpatialCoordinate(mesh)
n = FacetNormal(mesh)

# Define Taylor--Hood function space W
V = VectorFunctionSpace(mesh, "Lagrange" , 2)
Q = FunctionSpace(mesh , "Lagrange", 1)
P2 = VectorElement("Lagrange", triangle, 2)
P1 = FiniteElement("Lagrange", triangle, 1)
TH = MixedElement([P2, P1])
W = FunctionSpace(mesh, TH)

# Define Function and TestFunction(s)
w = Function(W)
(u, p) = split(w)
(v, q) = TestFunctions(W)

# Define bcs
p0 = (1.0 - x[0])  # or Expression("1.0-x[0]", degree=1)
bcs = DirichletBC(W.sub(0), (0.0, 0.0),
                  "on_boundary && !(near(x[0], 0.0) || near(x[0], 1.0))")

# Define initial viscosity guess
nu_guess = 0.2

# Define actual viscosity
def nu(u):
    return 0.5*pow(grad(u)**2, 1.0/(2*(4-1)))

# Define variational form for guess
epsilon = sym(grad(u))
F = (2*nu_guess*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx\
    + p0*dot(v,n)*ds

# Solve problem for guess
solve(F == 0, w, bcs)

# Plot initial guesses
print("||u_guess||_0^2 = ", assemble(u**2*dx))
plot(u, title="Velocity initial guess")
plot(p, title="Pressure initial guess")

# Define actual problem
F = (2*nu(u)*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx\
    + p0*dot(v,n)*ds

# Solve actual problem
solve(F == 0, w, bcs)

# Plot solutions
print("||u||_0^2 = ", assemble(u**2*dx))
plot(u, title="Velocity")
plot(p, title="Pressure")
interactive()
