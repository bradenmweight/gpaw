from ase import Atoms
from gpaw import GPAW, FermiDirac, Mixer
from gpaw import PW
from gpaw.hybrids.energy import non_self_consistent_energy as nsc_energy

from ase.parallel import paropen

# This calculation is too heavy to run as an exercise!!

myresults = paropen('si.atom.pbe_and_exx_energies.txt', 'a')

# Plane wave cutoff
pwcutoff = 400.0

# Do the isolated calculation

for L in [6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]:
    isolated_silicon = Atoms(
        symbols=['Si'],
        positions=[[0.5 * L, 0.5 * L, 0.5 * L]],
        cell=[L + 0.1, L, L - 0.1],
        pbc=(1, 1, 1))

    isolated_calc = GPAW(
        mode=PW(pwcutoff, force_complex_dtype=True),
        xc='PBE',
        txt='si_isolated_pbe.txt',
        occupations=FermiDirac(0.01),
        spinpol=True,
        hund=True,
        # convergence={'density': 1.e-6},
        mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0))

    isolated_silicon.calc = isolated_calc

    e0_isolated_pbe = isolated_silicon.get_potential_energy()
    isolated_calc.write('si.pbe+exx.isolated.gpw', mode='all')

    # Now the exact exchange
    si_isolated_exx = nsc_energy('si.pbe+exx.isolated.gpw',
                                 'EXX').sum()

    s = str(L)
    s += ' '
    s += str(e0_isolated_pbe)
    s += ' '
    s += str(si_isolated_exx)
    s += '\n'
    myresults.write(s)
