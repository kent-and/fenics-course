
from fenics import *
from fileplot import plot

def test_fileplot():
    mesh = UnitSquareMesh(10, 10)
    plot(mesh, title="mesh")
    V = FunctionSpace(mesh, "Lagrange", 1)
    V2 = FunctionSpace(mesh, "Lagrange", 2)
    W = VectorFunctionSpace(mesh, "Lagrange", 1)

    f = Function(V)
    f.interpolate(Expression("x[0]*x[0]*x[0]"))
    plot(f, title="f")

    u = Function(V2, name="u")
    u.interpolate(Expression("x[0]*x[0]*x[0]"))
    plot(u)

    g = f**2
    plot(g, title="g")

    plot(1.0, title="one", mesh=mesh)

    w = Function(W)
    w.interpolate(Expression(("sin(2.0*DOLFIN_PI*x[0])", "cos(2.0*DOLFIN_PI*x[1])")))
    plot(w, title='w')

    plot(grad(g), title='grad_g')

if __name__ == "__main__":
    test_fileplot()
