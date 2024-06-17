import sys

from ase import Atoms
from ase.build.connected import connected_indices
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
from ase.data.s22 import data
from ase.parallel import paropen

from gpaw import GPAW, FermiDirac
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii
from gpaw.utilities.adjust_cell import adjust_cell

try:
    from dftd4 import D4_model
except ModuleNotFoundError:
    pass

h = 0.18
box = 4.

xc = 'TS09'
if len(sys.argv) > 1:
    xc = sys.argv[1]

f = paropen('energies_' + xc + '.dat', 'w')
print('# h=', h, file=f)
print('# box=', box, file=f)
print('# molecule E[1]  E[2]  E[1+2]  E[1]+E[2]-E[1+2]', file=f)
for molecule in data:
    print(molecule, end=' ', file=f)
    ss = Atoms(data[molecule]['symbols'],
               data[molecule]['positions'])
    # split the structures
    s1 = ss[connected_indices(ss, 0)]
    s2 = ss[connected_indices(ss, -1)]
    assert len(ss) == len(s1) + len(s2)
    if xc == 'TS09' or xc == 'TPSS' or xc == 'M06-L' or xc == 'dftd4':
        c = GPAW(mode='fd', xc='PBE', h=h, nbands=-6,
                 occupations=FermiDirac(width=0.1))
    else:
        c = GPAW(mode='fd', xc=xc, h=h, nbands=-6,
                 occupations=FermiDirac(width=0.1))
    E = []
    for s in [s1, s2, ss]:
        s.calc = c
        adjust_cell(s, box, h=h)
        if xc == 'TS09':
            s.get_potential_energy()
            cc = vdWTkatchenko09prl(HirshfeldPartitioning(c),
                                    vdWradii(s.get_chemical_symbols(), 'PBE'))
            s.calc = cc
        elif xc == 'dftd4':
            s.get_potential_energy()
            cc = D4_model(xc='PBE', calc=c)
            s.calc = cc
        if xc == 'TPSS' or xc == 'M06-L':
            ene = s.get_potential_energy()
            ene += c.get_xc_difference(xc)
            E.append(ene)
        else:
            E.append(s.get_potential_energy())
    print(E[0], E[1], E[2], end=' ', file=f)
    print(E[0] + E[1] - E[2], file=f)
    f.flush()
f.close()
