import glob

import Numeric as num
from gpaw.utilities import equal
from ASE.ChemicalElements.symbol import symbols

from gpaw.setup import Setup
from gpaw.xc_functional import XCFunctional


##for xcname in ['LDA', 'PBE', 'revPBE']:
for xcname in ['LDA']:
    xcfunc = XCFunctional(xcname)
    for symbol in []:#'H']:#symbols:
        try:
            s = Setup(symbol, xcfunc, lmax=2)
        except IOError:
            continue
        print s.D_sp[0]
        e_kinetic = s.Kc + num.dot(s.D_sp[0], s.K_p)
        e_electrostatic = s.M + num.dot(s.D_sp[0], s.M_p) + \
                          num.dot(s.D_sp[0], num.dot(s.M_pp, s.D_sp[0]))
        H_sp = num.zeros(s.D_sp.shape, num.Float)
        e_xc = s.xc.calculate_energy_and_derivatives(s.D_sp, H_sp)
        print e_kinetic, e_electrostatic, e_xc
        equal(e_kinetic, 0.0, 2e-6)
        equal(e_electrostatic, 0.0, 0.011)
        equal(e_xc, 0.0, 5e-6)



