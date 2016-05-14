import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree
from gpaw.jellium import JelliumSlab
from gpaw.poisson import PoissonSolver
from gpaw.dipole_correction import DipoleCorrection
from gpaw import GPAW, Mixer
from gpaw.test import equal

rs = 5.0 * Bohr  # Wigner-Seitz radius
h = 0.24          # grid-spacing
a = 8 * h        # lattice constant
v = 3 * a        # vacuum
L = 10 * a       # thickness
k = 6           # number of k-points (k*k*1)

ps = PoissonSolver()
dc = DipoleCorrection(ps, 2)

ne = a**2 * L / (4 * np.pi / 3 * rs**3)

bc = JelliumSlab(ne, z1=v, z2=v + L)

surf = Atoms(pbc=(True, True, False),
             cell=(a, a, v + L + v))
surf.calc = GPAW(background_charge = bc,
                 poissonsolver = dc,
                 xc='LDA_X+LDA_C_WIGNER',
                 eigensolver='dav',
                 charge=-ne,
                 kpts=[k, k, 1],
                 h=h,
                 maxiter=300,
                 convergence={'density': 0.001},
                 mixer=Mixer(0.03, 7, 100),
                 nbands=int(ne / 2) + 15,
                 txt='surface.txt')
e = surf.get_potential_energy()

efermi = surf.calc.get_fermi_level()
# Get (x-y-averaged) electrostatic potential
# Must collect it from the CPUs
# https://listserv.fysik.dtu.dk/pipermail/gpaw-users/2014-January/002524.html
v = (surf.calc.hamiltonian.finegd.collect(surf.calc.hamiltonian.vHt_g,
                                     broadcast = True) * Hartree).mean(0).mean(0)
v = (surf.calc.hamiltonian.vHt_g * Hartree).mean(0).mean(0)

# Get the work function
phi1 = v[-1] - efermi

#print(phi1)
#equal(phi1, 2.70417600672, 1e-5)
# Reference value: Lang and Kohn, 1971, Theory of Metal Surfaces - Work function
# r_s = 5, work function = 2.73 eV

