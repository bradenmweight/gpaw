from ase.structure import molecule
from ase.parallel import parprint

from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.analyse.hirshfeld import HirshfeldDensity, HirshfeldPartitioning
from gpaw.analyse.wignerseitz import WignerSeitz
from gpaw.test import equal
from gpaw.mpi import rank

h = 0.3
gpwname = 'H2O' + str(h) + '.gpw'

def run(lastres=[]):
    results = []

    # Hirshfeld ----------------------------------------

    if 1:

        hd = HirshfeldDensity(calc)

        # check for the number of electrons
        expected = [[None, 10],
                    [[0, 1, 2], 10],
                    [[1, 2], 2],
                    [[0], 8],
                    ]
        #expected = [[[0, 2], 9], ]
        #expected = [[None, 10], ]
        for result in expected:
            indicees, result = result
            full, gd = hd.get_density(indicees)
            parprint('indicees', indicees, end=': ') 
            parprint('result, expected:', gd.integrate(full), result)
            equal(gd.integrate(full), result, 1.e-8)

        hp = HirshfeldPartitioning(calc)
        vr = hp.get_effective_volume_ratios()
        parprint('Hirshfeld:', vr)
        if len(lastres):
            equal(vr, lastres.pop(0), 1.e-10)
        results.append(vr)

    # Wigner-Seitz ----------------------------------------

    if 1:
        ws = WignerSeitz(calc.density.finegd, mol, calc)

        vr = ws.get_effective_volume_ratios()
        parprint('Wigner-Seitz:', vr)
        if len(lastres):
            equal(vr, lastres.pop(0), 1.e-10)
        results.append(vr)

    return results

# calculate
parprint('fresh:')
mol = Cluster(molecule('H2O'))
mol.minimal_box(3, h=h)
calc = GPAW(nbands=6,
            h = h, 
            txt=None)
calc.calculate(mol)
calc.write(gpwname)
lastres = run()

# load previous calculation
parprint('reloaded:')
calc = GPAW(gpwname, txt=None)
mol = calc.get_atoms()
run(lastres)
