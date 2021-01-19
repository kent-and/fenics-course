from fenics import *

# Create mesh
mesh  = UnitSquareMesh(30, 30)

# Create function space
Vele = VectorElement("Lagrange", triangle, 2)
Pele = FiniteElement("Lagrange", triangle, 1)
W = FunctionSpace(mesh, MixedElement([Vele, Pele]))

# Create boundary conditions
def noslip_boundary(x):
    return near(x[1], 0.0) or near(x[1], 1.0)

def inflow_boundary(x):
    return near(x[0], 0.0)

def outflow_boundary(x):
    return near(x[0], 1.0)

bcs = [DirichletBC(W.sub(0), (0, 0), noslip_boundary),
       DirichletBC(W.sub(1),  1,     inflow_boundary),
       DirichletBC(W.sub(1),  0,     outflow_boundary)]

# Create forms
f    = Constant((0, 0))
u, p = TrialFunctions(W)
v, q = TestFunctions(W)
a    = inner(grad(u), grad(v))*dx + dot(grad(p), v)*dx + div(u)*q*dx
L    = dot(f, v)*dx

# Compute solution
w = Function(W)
solve(a == L, w,bcs)

# Plot solution
u, p = w.split()
plot(u, title="u")
plot(p, title="p")
interactive()
