from fenics import *

# Create mesh
mesh  = UnitSquareMesh(30, 30)

# Choose degree of velocity-pressure spaces
#k, l = 1, 0
k, l = 1, 1
#k, l = 2, 1
#k, l = 3, 2

# Create function spaces
Vele = VectorElement("Lagrange", triangle, k)
Pele = FiniteElement("DG" if l == 0 else "Lagrange", triangle, l)
W = FunctionSpace(mesh, MixedElement([Vele, Pele]))

# Create boundary conditions
def noslip_boundary(x):
    return near(x[1], 0.0)

def slip_boundary(x):
    return near(x[1], 1.0)

bcs = [DirichletBC(W.sub(0), Constant((0, 0)), noslip_boundary),
       DirichletBC(W.sub(0), Constant((1, 0)), slip_boundary)]

# Create forms
f    = Constant((0, 0))
u, p = TrialFunctions(W)
v, q = TestFunctions(W)
R    = (inner(grad(u), grad(v)) - div(v)*p + div(u)*q - inner(f, v))*dx
a    = lhs(R)
L    = rhs(R)

# Compute solution
w = Function(W)
solve(a == L, w,bcs)

# Plot solution
u, p = w.split()
plot(u, title="u")
plot(p, title="p")

# Save solution to file
u_file = File("u_sol_P%dP%d.pvd" % (k,l))
u_file << u
p_file = File("p_sol_P%dP%d.pvd" % (k,l))
p_file << p

interactive()
