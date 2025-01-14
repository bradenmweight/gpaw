from ase import Atoms
from gpaw import GPAW
from gpaw import MixerSum

L = 6
name = 'N2'
a = Atoms('N2',
          [(L / 2 + 1.098 / 2, L / 2, L / 2),
           (L / 2 - 1.098 / 2, L / 2, L / 2)],
          cell=(L, L, L), pbc=False)

calc = GPAW(mode='fd',
            h=0.22,
            xc='PBE',
            convergence={'eigenstates': 1.0e-7},
            txt=name + '_PBE.txt',
            eigensolver='rmm-diis')
a.calc = calc
e_n2 = a.get_potential_energy()
n2t = calc.get_xc_difference('TPSS')
n2rt = calc.get_xc_difference('revTPSS')

a.calc = calc.new(xc='TPSS', txt=name + '_TPSS.txt')
e_n2t = a.get_potential_energy()

a.calc = calc.new(xc='revTPSS', txt=name + '_revTPSS.txt')
e_n2rt = a.get_potential_energy()

name = 'N'
b = Atoms('N', [(L / 2, L / 2, L / 2)], magmoms=[3],
          cell=(L, L, L), pbc=False)

calc = calc.new(mixer=MixerSum(0.02, 5, 100),
                hund=True,
                txt=name + '_PBE.txt')
b.calc = calc
e_n = b.get_potential_energy()
nt = calc.get_xc_difference('TPSS')
nrt = calc.get_xc_difference('revTPSS')

b.calc = calc.new(xc='TPSS', txt=name + '_TPSS.txt')
e_nt = b.get_potential_energy()

b.calc = calc.new(xc='revTPSS', txt=name + '_revTPSS.txt')
e_nrt = b.get_potential_energy()

print('Atm. Experiment  ', -228.5)
print('Atm. PBE         ', (e_n2 - 2 * e_n) * 23.06)
print('Atm. TPSS(nsc)   ', ((e_n2 + n2t) - 2 * (e_n + nt)) * 23.06)
print('Atm. TPSS        ', (e_n2t - 2 * e_nt) * 23.06)
print('Atm. revTPSS(nsc)', ((e_n2 + n2rt) - 2 * (e_n + nrt)) * 23.06)
print('Atm. revTPSS     ', (e_n2rt - 2 * e_nrt) * 23.06)
