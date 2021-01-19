from __future__ import print_function
from fenics import *
from ufl.classes import Expr

def project_to_cg1(x, mesh):
    # Create piecewise linear function space for this shape
    sh = x.shape()
    if sh == ():
        V = FunctionSpace(mesh, "Lagrange", 1)
    elif len(sh) == 1:
        V = VectorFunctionSpace(mesh, "Lagrange", 1, dim=sh[0])
    else:
        V = TensorFunctionSpace(mesh, "Lagrange", 1, shape=sh)
    return project(x, V)

def plot(x, title=None, mesh=None, **kwargs):
    """Dolfin plot() replacement writing CG1 function to
    a file called plot_*.pvd instead of showing visually.

    Import this in place of the dolfin plot to run programs
    in batch mode and still get the plot output for inspection.

    Use title argument to select a readable filename.

    Ignoring kwargs such as 'rescale', 'mode', 'interactive'.
    """

    if isinstance(x, Function):
        if title is None:
            title = x.name()
    if title is None:
        print("Missing title for plot, using 'untitled'.")
        title = 'untitled'
    filename = "plot_%s.pvd" % title

    if isinstance(x, (int,float)):
        x = as_ufl(x)

    if isinstance(x, Mesh):
        value = x

    elif isinstance(x, Expr):
        if mesh is None:
            mesh = x.domain().data()
        if not isinstance(mesh, Mesh):
            print("Cannot find mesh in expression")
            return

        if (isinstance(x, Function)
            and x.function_space().ufl_element().family() == "Lagrange"
            and x.function_space().ufl_element().degree() == 1):
            value = x
        else:
            print("Projecting to piecewise linears.")
            value = project_to_cg1(x, mesh)

    else:
        print("Invalid type for plotting: %s" % x.__class__)
        return

    print("Writing to file %s" % filename)
    f = File(filename)
    f << value
