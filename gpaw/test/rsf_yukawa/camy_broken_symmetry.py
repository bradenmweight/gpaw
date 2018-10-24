"""Document unfriendly behaviour for calculations with broken box sym."""
from ase import Atoms
from gpaw import GPAW, KohnShamConvergenceError
from gpaw.xc.hybrid import HybridXC
from gpaw.occupations import FermiDirac
from gpaw.mixer import MixerDif
from gpaw.eigensolvers import RMMDIIS

tio2 = Atoms('TiO2', [(0, 0, 0), (0.66, 0.66, 1.34), (0.66, 0.66, -1.34)])
tio2.center(vacuum=4)
tio2.translate([0.01, 0.02, 0.03])

c = {'energy': 0.001, 'eigenstates': 3, 'density': 3}

# For broken symmetry calculations MixerSum, MixerSum2 fail,
# Mixer sometimes fail (went to wrong direction)
# MixerDifs work (in this case)

tio2.calc = GPAW(txt='TiO2-CAMY-B3LYP-BS.txt', xc=HybridXC('CAMY_B3LYP'),
                 eigensolver=RMMDIIS(), maxiter=42, mixer=MixerDif(),
                 convergence=c,
                 occupations=FermiDirac(width=0.0, fixmagmom=True))
tio2.set_initial_magnetic_moments([2.0, -1.0, -1.0])
try:
    e_tio2 = tio2.get_potential_energy()
except KohnShamConvergenceError:
    pass
# dissoziation energy
print(tio2.calc.scf.converged)
assert tio2.calc.scf.converged
