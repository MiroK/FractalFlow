from dolfin import *


def hamiltonian(mesh, potential):
    '''Build the system for eigenvalue problem

    -\Delta u + chi*u = lmbda*u in Omega
                    du/dn = 0 on boundary
    '''
    # The problem
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Potential
    chi = potential(V)
    # Save the potential
    File('chi.pvd') << chi

    # bc = DirichletBC(V, Constant(0), 'on_boundary')
    a = inner(grad(u), grad(v))*dx + inner(chi*u, v)*dx
    b = inner(u, v)*dx
    L = inner(Constant(0), v)*dx

    A, _ = assemble_system(a, L)
    B, _ = assemble_system(b, L)

    # Return for slepc4py
    A, B = map(lambda x: as_backend_type(x).mat(), (A, B))

    return A, B, V

            
# --------------------------------------------------------------------

if __name__ == '__main__':
    from distance import dof_chi
    from solver import eigensolve
    
    mesh_file = 'levy.xml'
    curve_file = 'levy_facet_region.xml'

    mesh = Mesh(mesh_file)
    f = MeshFunction('size_t', mesh, curve_file)
    # Potential computed once the space is given (so is a function of V)
    potential = lambda space, f=f: dof_chi(f, space)

    A, B, V = hamiltonian(mesh, potential)
    # Compute eigenpairs
    pairs = eigensolve(A, B, V, params={})

    # Save
    for i, (w, v) in enumerate(pairs):
        print '%d-th smallest eigenvalue is %g' % (i, w)
        File('vector_%d.pvd' % i) << v

