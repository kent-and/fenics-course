from __future__ import print_function
"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = exp(-100(x^2 + y^2))

and homogeneous Dirichlet boundary conditions.

Note that we use a simplified error indicator, ignoring
edge (jump) terms and the size of the interpolation constant.
"""

# Copyright (C) 2008 Rolv Erlend Bredesen
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Anders Logg 2008-2011
#
# First added:  2008-04-03
# Last changed: 2013-08-16

from fenics import *
from numpy import array, sqrt
from math import pow

TOL = 5e-4           # Error tolerance
REFINE_RATIO = 0.50  # Refine 50 % of the cells in each iteration
MAX_ITER = 30        # Maximal number of iterations

# Plot convergence of norm
norms = []

# Create initial mesh
mesh = UnitSquareMesh(4, 4)
comm = mesh.mpi_comm()
source_str = "exp(-100.0*(pow(x[0], 2) + pow(x[1], 2)))"
source = eval("lambda x: " + source_str)
# Adaptive algorithm
for level in xrange(MAX_ITER):

    # Define variational problem
    V = FunctionSpace(mesh, "Lagrange", 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Expression(source_str, degree=4)
    a = dot(grad(v), grad(u))*dx
    L = v*f*dx

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, DomainBoundary())

    # Compute solution
    u = Function(V)
    A = assemble(a)
    b = assemble(L)
    bc.apply(A, b)
    solver = KrylovSolver("cg", "amg")
    solver.parameters["absolute_tolerance"] = 1e-15
    solver.parameters["relative_tolerance"] = 1e-15
    solver.solve(A, u.vector(), b)

    # Compute error indicators
    h = array([c.diameter() for c in cells(mesh)])
    K = array([c.volume() for c in cells(mesh)])
    R = array([abs(source([c.midpoint().x(), c.midpoint().y()])) for c in cells(mesh)])
    gamma = h*R*sqrt(K)

    # Compute error estimate
    E = sum([g*g for g in gamma])
    E = sqrt(MPI.sum(comm, E))
    print("Level %d: E = %g (TOL = %g)" % (level, E, TOL))

    # Check convergence
    #if E < TOL:
    #    info("Success, solution converged after %d iterations" % level)
    #    break

    # Mark cells for refinement
    cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
    gamma_0 = sorted(gamma, reverse=True)[int(len(gamma)*REFINE_RATIO)]
    gamma_0 = MPI.max(comm, gamma_0)
    for c in cells(mesh):
        cell_markers[c] = gamma[c.index()] > gamma_0

    # Refine mesh
    mesh = refine(mesh, cell_markers)

    n = norm(u)
    print("Norm: %.16g\n" % n)
    norms.append(n)

    # Plot mesh
    #plot(mesh)

from pylab import *
plot(norms, "-o")

print("")
for n in norms:
    print("%.16g" % n)

show()
