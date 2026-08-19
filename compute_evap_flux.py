"""
Uses the finite element method (implemented in legacy FEniCS) to
compute the evaporation flux.  The code solves for the dimensionless
vapour concentration field by solving the 3D Laplace equation.

The length of the beam has been non-dimensionalised to one.  The
non-dim width and height of the beam are W / L and h_s / L.  The
coordinates of the wetted surface of the substrate are given by
-1/2 < x < 1/2, -1/2 < y 1/2, and z = 0. 

The beam is contained in a box of dimension H x H x H.  Due to symmetry,
the equations are solved in the region 0 < x < H / 2, 0 < y < H / 2,
and -H < z < H

"""
from fenics import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

"""
Define the output files for the concentration u and
the flux q
"""

file = File('fem_soln/u.pvd')
file_q = File('fem_soln/q.pvd')


# non-dimensional boundary of the box
H = 10

"""
Import mesh
"""
mesh_name = 'meshes/full'
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), mesh_name+'.h5', "r")
hdf.read(mesh, "/mesh", False)
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")
bdry = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
hdf.read(bdry, "/boundaries")

# id of the wetted surface/film
film = 999

# re-define the measures
dx = Measure('dx', mesh, subdomain_data=subdomains)
print(assemble(Constant(1) * dx))

"""
Function spaces and functions
"""

# function spaces
P2 = FiniteElement('CG', mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, P2)

# functions
u = Function(V)
u_ = TestFunction(V)

"""
Boundary conditions.  Impose c = 1 at the wetted surface and
c = 0 in the far field.  Symmetry conditions automatically imposed
in the weak form.
"""

# define the far-field boundaries
def sides(x):
    return near(abs(x[2]), H/2) or near(x[1], H/2) or near(x[0], H/2)


# Dirichlet boundary conditions
BC = [
    DirichletBC(V, Constant(1), bdry, film),
    DirichletBC(V, Constant(0), sides)
    ]


"""
Build the weak form and its Jacobian
"""
F = -inner(grad(u), grad(u_)) * dx
J = derivative(F, u)


"""
Set up solver
"""
problem = NonlinearVariationalProblem(F, u, BC, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['nonlinear_solver'] = 'snes'
prm['snes_solver']['linear_solver'] = 'gmres'
prm['snes_solver']['absolute_tolerance'] = 1E-8
prm['snes_solver']['relative_tolerance'] = 1E-7
prm['snes_solver']['maximum_iterations'] = 25


"""
Solve the problem
"""
solver.solve()

file << u


"""
Compute flux by projecting onto DG0 elements
"""

W = VectorFunctionSpace(mesh, "DG", 0)
g = TrialFunction(W)
v = TestFunction(W)

a = inner(g, v)*dx
L = inner(grad(u), v)*dx

A = assemble(a)
b = assemble(L)

q = Function(W)

q.rename('q', 'q')


solver = PETScKrylovSolver("cg", "hypre_amg")
solver.solve(A, q.vector(), b)

file_q << q