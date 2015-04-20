from ase.dft.kpoints import ibz_points, get_bandpath
from ase.parallel import paropen
from gpaw import GPAW
import numpy as np

layer = GPAW('Fe_gs.gpw', txt=None).atoms

points = ibz_points['bcc']
G = points['Gamma']
H = points['H']
P = points['P']
N = points['N']
H_z = [H[0], -H[1], -H[2]]
G_yz = [2*H[0], 0.0, 0.0]

kpts, x, X = get_bandpath([G, H, G_yz], layer.cell, npoints=500)
calc = GPAW('Fe_gs.gpw',
            kpts=kpts,
            symmetry='off',
            txt='Fe_bands.txt',
            parallel={'band': 1})
calc.diagonalize_full_hamiltonian()

calc.write('Fe_bands.gpw', mode='all')

f = paropen('Fe_kpath.dat', 'w')
for k in x:
    print >> f, k
f.close()

f = paropen('Fe_highsym.dat', 'w')
for k in X:
    print >> f, k
f.close()