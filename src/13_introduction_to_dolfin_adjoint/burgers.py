from fenics import *

n = 30
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, "Lagrange", 2)

ic = Expression(("sin(2*pi*x[0])", "cos(2*pi*x[1])"), degree=4)
u = project(ic, V)
u_next = Function(V)
v = TestFunction(V)

nu = Constant(0.0001)

timestep = Constant(0.01)

F = (inner((u_next - u)/timestep, v)
     + inner(grad(u_next)*u_next, v)
     + nu*inner(grad(u_next), grad(v)))*dx

bc = DirichletBC(V, (0.0, 0.0), "on_boundary")

t = 0.0
end = 0.1
while (t <= end):
    solve(F == 0, u_next, bc)
    plot(u_next)
    u.assign(u_next)
    t += float(timestep)

interactive()
