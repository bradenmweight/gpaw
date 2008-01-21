import sys
class Out:
    def write(self, x):
        sys.__stdout__.write(x)
        raise RuntimeError('not silent')

out, err = sys.stdout, sys.stderr
sys.stdout = sys.stderr = Out()

try:
    from gpaw import Calculator
    from ase import *

    a = 5.0
    h = 0.2
    calc = Calculator(h=h, nbands=1, kpts=(1, 1, 1),
                      width=1e-9,
                      xc='PBE',
                      txt=None)
    hydrogen = Atoms([Atom('H', (a / 2, a / 2, a / 2), magmom=0)],
                     cell=(a, a, a),
                     calculator=calc)
    f = hydrogen.get_forces()
except:
    sys.stdout = out
    sys.stderr = err
    raise

sys.stdout = out
sys.stderr = err
