__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-19"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

# Modified by Marie E. Rognes

# Last changed: 2012-06-25

from fenics import *

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

mesh = UnitSquareMesh(6, 4)
V = FunctionSpace(mesh, "Lagrange", 1)

u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
boundary = "on_boundary"  # Alt 1
#boundary = Boundary()    # Alt 2
bc = DirichletBC(V, u0, boundary)

f = Constant(-6.0)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

u = Function(V)
solve(a == L, u, bc)

plot(u)
interactive()
