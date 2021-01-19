# This is a simple benchmark problem that can be used to test the
# speed of students laptops and FEniCS installations. The size of
# the linear system is 2,146,689 x 2,146,689.
#
# For reference, this benchmark takes 18.3 seconds on my laptop.
#
# Anders Logg 2013-08-14

from fenics import *

set_log_level(PROGRESS)

# Size of problem
N = 128

# Check which preconditioners are available
if ("amg", "Algebraic multigrid") in krylov_solver_preconditioners():
    pc = "amg"
else:
    print("Warning: No AMG preconditioner available, using ILU")
    pc = "ilu"

# Set up problem
mesh = UnitCubeMesh(N, N, N)
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = u*v*dx + dot(grad(u), grad(v))*dx
L = Constant(100.0)*v*dx

# Time finite element assembly
tic()
A = assemble(a)
b = assemble(L)
t1 = toc()

# Time linear solve (iterative solver)
tic()
x = Vector()
solve(A, x, b, "cg", pc)
t2 = toc()

# Report timings
print("")
print("Finite element assembly:  %.3g s" % t1)
print("Linear solve (iterative): %.3g s" % t2)
print("")
print("Total time:               %.3g s" % (t1 + t2))
