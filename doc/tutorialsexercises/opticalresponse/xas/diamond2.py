from gpaw import GPAW
from gpaw.xas import RecursionMethod

name = 'diamond333_hch'

calc = GPAW(name + '.gpw',
            kpts=(6, 6, 6),
            txt=name + '_rec.txt')
calc.initialize()
calc.set_positions()

r = RecursionMethod(calc)
r.run(600)

r.run(1400,
      inverse_overlap='approximate')

r.write(name + '_600_1400a.rec',
        mode='all')
