from ase import Atoms
from ase.build.connected import connected_indices
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
from ase.data.s22 import data
from ase.parallel import parprint

from gpaw import GPAW, FermiDirac
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii
from gpaw.utilities.adjust_cell import adjust_cell

h = 0.25
box = 3.0

molecule = 'Adenine-thymine_complex_stack'

Energy = {
    'PBE': [],
    'vdW-DF': [],
    'TS09': []}

for molecule in ['Adenine-thymine_complex_stack']:
    ss = Atoms(data[molecule]['symbols'],
               data[molecule]['positions'])

    # split the structures
    s1 = ss[connected_indices(ss, 0)]
    s2 = ss[connected_indices(ss, -1)]
    assert len(ss) == len(s1) + len(s2)
    calc_params = dict(mode='fd', h=h, nbands=-6,
                       occupations=FermiDirac(width=0.1), txt=None)
    c = GPAW(**calc_params, xc='PBE')
    cdf = GPAW(**calc_params, xc='vdW-DF')

    for s in [s1, s2, ss]:
        s.calc = c
        adjust_cell(s, box, h=h)
        Energy['PBE'].append(s.get_potential_energy())
        cc = vdWTkatchenko09prl(HirshfeldPartitioning(c),
                                vdWradii(s.get_chemical_symbols(), 'PBE'))
        s.calc = cc
        Energy['TS09'].append(s.get_potential_energy())

        s.calc = cdf
        Energy['vdW-DF'].append(s.get_potential_energy())

    parprint('Coupled cluster binding energy',
             -data[molecule]['interaction energy CC'] * 1000, 'meV')
    for xc in ['PBE', 'vdW-DF', 'TS09']:
        ene = Energy[xc]
#        print xc, 'energies', ene
        parprint(xc, 'binding energy',
                 (ene[0] + ene[1] - ene[2]) * 1000, 'meV')
