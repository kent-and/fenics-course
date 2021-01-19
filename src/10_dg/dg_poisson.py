from fenics import *

n = 32
mesh = UnitSquareMesh(n, n)
n = FacetNormal(mesh)
h = avg(mesh.hmin())

V = FunctionSpace(mesh, "DG", 0)
u = TrialFunction(V)
v = TestFunction(V)

alpha = Constant(10.0)
f = Constant(1.0)

a_K = inner(grad(u), grad(v))*dx
a_int = (- inner(avg(grad(u)), jump(v, n)) - inner(jump(u, n), avg(grad(v)))
         + alpha/h*dot(jump(u, n), jump(v, n)))*dS
a_ext = (- inner(dot(grad(u), n), v) - inner(dot(grad(v), n), u)
         + alpha/h*u*v)*ds

a = a_K + a_int + a_ext
L = f*v*dx

u = Function(V)
solve(a == L, u)
plot(u, interactive=True)
