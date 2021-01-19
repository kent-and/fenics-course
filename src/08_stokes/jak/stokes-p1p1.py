from fenics import *

# Read in mesh
mesh = Mesh("dolfin-2.xml.gz")

# Create mixed function space
Vele = VectorElement("Lagrange", triangle, 1)
Pele = FiniteElement("Lagrange", triangle, 1)
W = FunctionSpace(mesh, MixedElement([Vele, Pele]))

# Define boudary domains
def inflow_boundary(x):
    return near(x[0],1.0)

def outflow_boundary(x):
    return near(x[0],0.0)

def noslip_boundary(x, on_boundary):
    tol = 1.0e-14
    return (x[1] < tol or x[1] > 1 - tol
           or (on_boundary and not inflow_boundary(x) and not outflow_boundary(x)))

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

# Stabilization parameters
#beta = Constant(0.0)
beta = Constant(0.1)
#beta = Constant(100.0)
h = CellSize(mesh)
delta = beta*h*h

# Define variational form
a = (inner(grad(u), grad(v)) - div(v)*p - q*div(u)
    - delta*inner(grad(p), grad(q)))*dx
L = inner(f, v)*dx - delta*inner(f, grad(q))*dx

# Compute solution
w = Function(W)
solve(a == L, w, bcs)

# Split and plot solution
(u, p) = w.split()
plot(u)
plot(p)

# Save solution in VTK format
File("velocity.pvd") << u
File("pressure.pvd") << p

interactive()
