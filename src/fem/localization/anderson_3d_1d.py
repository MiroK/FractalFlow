from dolfin import *
from distance import dof_chi
from solver import eigensolve
from hamiltonian import hamiltonian
import os

which_curve = 'line'

if which_curve == 'anna':
    mesh = UnitCubeMesh(32, 32, 32)
    f = MeshFunction('size_t', mesh, 1, 0)
    CompiledSubDomain('near(x[0], 0.5) && near(x[1], 0.5)').mark(f, 1)

else:
    dir = '/mn/sarpanitu/eksterne-u1/mirok/Documents/Programming/enumath3d1d/src/mesh'
    
    mesh_file = os.path.join(dir, 'vasculature_volmax2000_mesh.xml')
    curve_file = os.path.join(dir, 'vasculature_volmax2000_markers.xml')

    mesh = Mesh(mesh_file)
    f = MeshFunction('size_t', mesh, curve_file)

# Potential computed once the space is given (so is a function of V)
potential = lambda space, f=f: dof_chi(f, space)

A, B, V = hamiltonian(mesh, potential)
info('dim(V) is %d' % V.dim())
# Compute eigenpairs
pairs = eigensolve(A, B, V, params={'-eps_nev': 10}, small_enough=-1)

# Save
for i, (w, v) in enumerate(pairs):
    info('%d-th smallest eigenvalue is %g' % (i, w))
    v.rename('f', '0')
    File('%s_vector_%d.pvd' % (which_curve, i)) << v

phi = potential(V)
phi.rename('f', '0')
File('%s_potential.pvd' % which_curve) << phi

File('%s_curve.pvd' % which_curve) << f
