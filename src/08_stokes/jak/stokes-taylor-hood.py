from fenics import *

# Read in mesh
mesh = Mesh("dolfin-2.xml.gz")

# Create mixed function space
Vele = VectorElement("Lagrange", triangle, 2)
Pele = FiniteElement("Lagrange", triangle, 1)
W = FunctionSpace(mesh, MixedElement([Vele, Pele]))

# Define boudary domains
def outflow_boundary(x):
    return near(x[0],0.0)

def inflow_boundary(x):
    return near(x[0],1.0) #or near(x[0],0)

def noslip_boundary(x, on_boundary):
    return on_boundary and not inflow_boundary(x) and not outflow_boundary(x)

# Define boundary conditions
inflow = Expression(("-sin(pi*x[1])", "0"), degree=4)
noslip = Constant((0, 0))
zero = Constant(0)

bc1 = DirichletBC(W.sub(0), inflow, inflow_boundary)
bc2 = DirichletBC(W.sub(0), noslip, noslip_boundary)
bc3 = DirichletBC(W.sub(1), zero, outflow_boundary)

# Collect boundary conditions
bcs = [bc1, bc2, bc3]

# Define test and trial functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Source term
f = Constant((0,0))

# Define variational form
a = (inner(grad(u), grad(v)) - div(v)*p - q*div(u))*dx
L = inner(f, v)*dx

# Compute solution
w = Function(W)
solve(a == L, w, bcs)

# Split and plot solution
(u,p) = w.split()
plot(u)
plot(p)

# Save solution in VTK format
File("velocity.pvd") << u
File("pressure.pvd") << p

interactive()
