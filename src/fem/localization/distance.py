# Getting distance to curve
from dolfin import Function, FunctionSpace, SubsetIterator, as_backend_type


def dof_chi(f, V, marker=1):
    '''
    Let f be a an edge function. Build a function in V which is such 
    that all dofs where f == marker are 1 and are 0 otherwise.
    '''
    tdim = f.mesh().topology().dim()
    assert f.dim() == 1

    assert V.mesh().id() == f.mesh().id()
    assert V.ufl_element().family() == 'Lagrange'

    mesh = V.mesh()
    mesh.init(1, tdim)
    mesh.init(tdim, 1)
    e2c = mesh.topology()(1, tdim)
    c2e = mesh.topology()(tdim, 1)
    
    dm = V.dofmap()
    first, last = dm.ownership_range()
    n = last - first

    is_local = lambda i, f=first, l=last: f <= i+f < l
    
    # Triangle has 3 facets/edges, tet has 6 edges
    nedges = mesh.ufl_cell().num_edges()
    edge_dofs = [dm.tabulate_entity_closure_dofs(1, i) for i in range(nedges)]
    
    chi = Function(V)
    values = chi.vector().get_local()
    for edge in SubsetIterator(f, marker):
        edgei = edge.index()

        c = e2c(edgei)[0]  # Any cell
        dofs = dm.cell_dofs(c)

        local_edge = c2e(c).tolist().index(edgei)
        dofs = filter(is_local, dofs[edge_dofs[local_edge]])
        values[dofs] = 1.
    # Sync
    chi.vector().set_local(values)
    as_backend_type(chi.vector()).update_ghost_values()

    return chi


def cell_chi(f, marker=1):
    '''
    Let f be a facet function. Build a DG0 function which is such 
    that cell connected to f == marker are 1 and are 0 otherwise.
    '''
    assert f.mesh().topology().dim() == 2
    assert f.dim() == 1

    mesh = f.mesh()
    mesh.init(1, 2)
    f2c = mesh.topology()(1, 2)

    V = FunctionSpace(mesh, 'DG', 0)

    dm = V.dofmap()
    first, last = dm.ownership_range()
    n = last - first

    is_local = lambda i, f=first, l=last: f <= i+f < l
    
    chi = Function(V)
    values = chi.vector().get_local()
    for facet in SubsetIterator(f, marker):
        fi = facet.index()

        for c in f2c(fi):
            dofs = filter(is_local, dm.cell_dofs(c))
            values[dofs] = 1.
    # Sync
    chi.vector().set_local(values)
    as_backend_type(chi.vector()).update_ghost_values()

    return chi


def distance_chi(curve_f, V, p, eps, solver_params):
    '''
    Make a P1 function f which is a distance function from curve_f.
    That is 
      
      f(x) = min_{y in mesh}dist(x, y, p)

    where
    
      dist(x, y) = ||x-y||_p = (|x_0 - y_0|^{p} + |x_1 - y_1|^{p})^{1/p}
    '''
    # NOTE: this might require regulatization -> eps
    #       obviously diferent p will change the metric
    pass

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    
    mesh = UnitSquareMesh(32, 32)

    curve_f = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    CompiledSubDomain('near(x[0], x[1])').mark(curve_f, 1)

    V = FunctionSpace(mesh, 'CG', 1)

    # f = dof_chi(curve_f, V)
    f = cell_chi(curve_f)

    out = XDMFFile('f.xdmf')
    out.write(f)
