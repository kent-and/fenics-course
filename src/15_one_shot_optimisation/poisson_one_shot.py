from fenics import *

# Define discrete Functionspace
mesh = UnitSquareMesh(100, 100)
U = FiniteElement("Lagrange", triangle, 1) # Finite element for forward PDE space
V = FiniteElement("Lagrange", triangle, 1) # Finite element for adjoint PDE space
M = FiniteElement("DG", triangle, 0) # Finite lement for control space
W = FunctionSpace(mesh, MixedElement([U, V, M]))

# Define Functions
w = Function(W)
u, lmbd, f = split(w) # Solution functions

x = TestFunction(W)  # Test functions


# Define variational problem
a = inner(grad(u), grad(lmbd))*dx
L = f*lmbd*dx

# Define functional
ud = Expression("sin(pi*x[0])*sin(pi*x[1])", degree=4)   # Desired temperature profile
alpha = Constant(1e-6) # Regularization parameter
J = (u-ud)**2*dx + alpha*f**2*dx

# Define boundary conditions
bc_u    = DirichletBC(W.sub(0), 0.0, "on_boundary")
bc_lmbd = DirichletBC(W.sub(1), 0.0, "on_boundary")
bcs = [bc_u, bc_lmbd]

# Derive optimality conditions
lagrang = J + a + L
kkt = derivative(lagrang, w, x)


# Solve Poisson problem
solve(kkt == 0, w, bcs)

plot(w[0], title="Temperature")
plot(w[2], title="Control")
plot(ud-w[0], title="Temperature difference")
interactive()
