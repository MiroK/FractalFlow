from petsc4py import PETSc
from slepc4py import SLEPc


def eigensolve(A, B, V, params, small_enough=5000):
    '''Solve A*u = lB*u returning the eigenpairs'''
    # Do small enought exacty
    if V.dim() < small_enough: return exact_eigensolve(A, B, V, params)
    
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
