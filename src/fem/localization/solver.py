from scipy.linalg import eigh
from scipy.sparse import csr_matrix
import numpy as np
from dolfin import Function, PETScVector
from petsc4py import PETSc
from slepc4py import SLEPc


def exact_eigensolve(A, B, V, params):
    '''A direct solver intended to run in serial'''
    assert A.comm.size == 1

    A = csr_matrix(A.getValuesCSR()[::-1], shape=A.size)
    B = csr_matrix(B.getValuesCSR()[::-1], shape=B.size)
    
    eigw, eigv = eigh(A.todense(), B.todense())
    sort_idx = np.argsort(eigw)

    # Fall back to 10 eigenpair
    nev = params.get('-eps_type', 10)
    eigw = eigw[sort_idx[:nev]]
    eigv = (eigv.T)[sort_idx[:nev]]

    eigenpairs = []
    for w, v in zip(eigw, eigv):
        f = Function(V)
        f.vector().set_local(v)

        eigenpairs.append((w, f))
    return eigenpairs


def eigensolve(A, B, V, params, small_enough=5000):
    '''Solve A*u = lB*u returning the eigenpairs'''
    # Do small enought exacty
    if V.dim() < small_enough:
        print 'Using scipy as dim(V) is %d' % V.dim()
        return exact_eigensolve(A, B, V, params)
    
    # NOTE: you configure this from command line
    # Here are some defaults
    my_params = {'-eps_tol': 1E-6,         # cvrg tolerance
                 '-eps_max_it': 10000,      
                 '-eps_smallest_magnitude': 'none',  # which eigenvalues
                 '-eps_nev': 3,                      # How many
                 '-eps_monitor': 'none',
                 '-eps_type': 'krylovschur'}
    
    for key, value in my_params.items():
        if key not in params:
            params[key] = value

    opts = PETSc.Options()
    for key, value in params.items():
        opts.setValue(key, None if value == 'none' else value)

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setOperators(A ,B)
    # type is -eps_type
    E.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    E.setFromOptions()

    # FIXME: spectral transform and precond to accelerate things
    E.solve()

    its = E.getIterationNumber()
    nconv = E.getConverged()
    assert nconv > 0

    eigenpairs = []
    for i in range(nconv):
        eigv = A.createVecLeft()
        eigw = E.getEigenpair(i, eigv).real

        eigenpairs.append((eigw, Function(V, PETScVector(eigv))))
    return eigenpairs
