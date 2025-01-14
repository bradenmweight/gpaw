from ase.build import mx2

from gpaw import GPAW
from gpaw.unfold import Unfold, find_K_from_k

a = 3.184
PC = mx2(a=a).get_cell(complete=True)
bp = PC.get_bravais_lattice().bandpath('MKG', npoints=48)
x, X, _ = bp.get_linear_kpoint_axis()

M = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]

Kpts = []
for k in bp.kpts:
    K = find_K_from_k(k, M)[0]
    Kpts.append(K)

calc_bands = GPAW('gs_3x3_defect.gpw').fixed_density(
    kpts=Kpts,
    symmetry='off',
    nbands=220,
    convergence={'bands': 200})

calc_bands.write('bands_3x3_defect.gpw', 'all')

unfold = Unfold(name='3x3_defect',
                calc='bands_3x3_defect.gpw',
                M=M,
                spinorbit=False)

unfold.spectral_function(kpts=bp.kpts, x=x, X=X,
                         points_name=['M', 'K', 'G'])
