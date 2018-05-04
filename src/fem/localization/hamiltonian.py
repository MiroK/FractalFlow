from slepc4py import SLEPc
from petsc4py import PETSc
from scipy.linalg import eigh
from dolfin import *



def hamiltonian(mesh, curve_f):
    '''Build the system for eigenvalue problem

    -\Delta u + chi*u = lmbda*u in Omega
                    du/dn = 0 on boundary
    '''
    # The problem
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Potential
    chi = curve_chi(curve_f, V)

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

    
    mesh = Mesh('../test.xml')

    f = MeshFunction('size_t', mesh, '../test_facet_region.xml')

    A, B, V = hamiltonian(mesh, f)
    E = eigensolver(A, B, V, params={})
    

    #V = FunctionSpace(mesh, 'CG', 1)

    #g = curve_chi(f, V)

    #f = XDMFFile('foo.xdmf')
    #f.write(g)

    # # Load mesh and the marking function
    # comm = mpi_comm_world()
    # h5 = HDF5File(comm, mesh_path, 'r')
    # mesh = Mesh()
    # h5.read(mesh, 'mesh', False)

    # curve_f = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    # h5.read(curve_f, 'facet')

