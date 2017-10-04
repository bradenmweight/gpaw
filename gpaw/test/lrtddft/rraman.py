from ase import Atoms, Atom
from ase.vibrations.placzek import Placzek

from gpaw import GPAW
from gpaw.lrtddft.kssingle import KSSingles
from gpaw.test import equal

txt = '-'
txt = None
load = True
load = False
xc = 'LDA'

R = 0.7  # approx. experimental bond length
a = 4.0
c = 5.0
H2 = Atoms([Atom('H', (a / 2, a / 2, (c - R) / 2)),
            Atom('H', (a / 2, a / 2, (c + R) / 2))],
           cell=(a, a, c))

calc = GPAW(xc=xc, nbands=3, spinpol=False, eigensolver='rmm-diis', txt=txt)
H2.set_calculator(calc)
H2.get_potential_energy()

gsname = exname = 'rraman'
pz = Placzek(H2, KSSingles, gsname=gsname, exname=exname,
             exkwargs={'eps':0.0, 'jend':1}
)
pz.run()

# check size
kss = KSSingles('rraman-d0.010.eq.ex.gz')
assert(len(kss) == 1)

pz = Placzek(H2, KSSingles, gsname=gsname, exname=exname, 
                   verbose=True,)
ai = pz.absolute_intensity(omega=5)[-1]
equal(ai, 299.955946107, 1e-3) # earlier obtained value
i = pz.intensity(omega=5)[-1]
equal(i, 7.82174195159e-05, 1e-11) # earlier obtained value
pz.summary(omega=5, method='frederiksen')
