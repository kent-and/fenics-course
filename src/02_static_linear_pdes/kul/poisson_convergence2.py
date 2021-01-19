from __future__ import print_function
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-19"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

# Modified by Martin Alnes
# Last changed: 2016-04-10

from fenics import *
import pylab as p
from plotslopes import *

def solve_poisson(q, n):

    # Create mesh and define function space
    mesh = UnitSquareMesh(n, n)
    x = SpatialCoordinate(mesh)
    V = FunctionSpace(mesh, "Lagrange", q)

    # Print size of system
    N = V.dim()
    print("n = %d num_dofs = %d" % (n, N))

    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x, on_boundary):
        return on_boundary

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)

    # Define rhs using Expression or symbolic expression
    f = Expression("2.0*DOLFIN_PI*DOLFIN_PI*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])", degree=4)
    #f = (2.0*pi**2) * sin(pi*x[0]) * sin(pi*x[1])

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx(degree=4)  # Note specification of quadrature degree for integration

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    # Compute error using symbolic exact solution
    u_exact = sin(pi*x[0]) * sin(pi*x[1])
    degree = 2*(q + 2)
    error = sqrt( assemble((u_exact - u)**2*dx(degree=degree)) )

    # Compute error using Expression and errornorm()
    #u_exact = Expression("sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])", degree=8)
    #error = errornorm(u_exact, u, degree_rise=2)

    return 1.0 / n, N, error

# Check convergence
h1, N1, e1 = zip(*[solve_poisson(1, n) for n in [2, 4, 8, 16, 32, 64, 128, 256]])
h2, N2, e2 = zip(*[solve_poisson(2, n) for n in [2, 4, 8, 16, 32, 64, 128, 256]])

def fmt(values):
    return "  ".join("%.4e" % v for v in values)
print("Convergence:")
print("h1:", fmt(h1))
print("h2:", fmt(h2))
print("N1:", fmt(N1))
print("N2:", fmt(N2))
print("e1:", fmt(e1))
print("e2:", fmt(e2))

# Plot results
enable_plot = False
if enable_plot:
    p.hold(True)
    p.loglog(h1, e1, 'g-o')
    p.loglog(h2, e2, 'b-o')
    slope_marker((h1[5], e1[5]), (2, 1))
    slope_marker((h2[5], e2[5]), (3, 1))
    p.legend(["P1", "P2"], loc="lower right")
    p.xlabel("h")
    p.ylabel("e")
    p.grid(True)
    p.show()
