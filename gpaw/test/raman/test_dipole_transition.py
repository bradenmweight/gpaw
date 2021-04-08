import numpy as np
import pytest

from ase.units import Bohr
from gpaw import GPAW
from gpaw.raman.dipoletransition import get_dipole_transitions
from gpaw.utilities.dipole import dipole_matrix_elements_from_calc
from gpaw.lrtddft.kssingle import KSSingles


def test_dipole_transition(gpw_files, tmp_path_factory):
    """Check dipole matrix-elements for H20."""
    calc = GPAW(gpw_files['h2o_lcao_wfs'])
    dip_svknm = get_dipole_transitions(calc.atoms, calc, savetofile=False).real
    assert dip_svknm.shape == (1, 3, 1, 6, 6)
    dip_vnm = dip_svknm[0, :, 0]

    for i in range(3):
        # d2 = np.real(dip_svknm[0, i, 0] * dip_svknm[0, i, 0].conj())
        assert(np.allclose(dip_svknm[0, i, 0] + dip_svknm[0, i, 0].T, 0.,
                           atol=1e-4))
        # assert(np.allclose(d2, d2.T, rtol=1e-3, atol=1e-5))

    # Check numerical value of a few elements - signs might change!
    assert 0.3265 == pytest.approx(abs(dip_svknm[0, 0, 0, 0, 3]), abs=1e-4)
    assert 0.1411 == pytest.approx(abs(dip_svknm[0, 0, 0, 2, 3]), abs=1e-4)
    assert 0.0987 == pytest.approx(abs(dip_svknm[0, 0, 0, 3, 4]), abs=1e-4)
    assert 0.3265 == pytest.approx(abs(dip_svknm[0, 0, 0, 3, 0]), abs=1e-4)
    assert 0.3889 == pytest.approx(abs(dip_svknm[0, 1, 0, 0, 1]), abs=1e-4)
    assert 0.3669 == pytest.approx(abs(dip_svknm[0, 2, 0, 0, 2]), abs=1e-4)

    # some printout for manual inspection, if wanted
    f = 6 * "{:+.4f} "
    for c in range(3):
        for i in range(6):
            print(f.format(*dip_vnm[c, i]))
        print("")

    # ------------------------------------------------------------------------
    # compare to utilities implementation
    uref = dipole_matrix_elements_from_calc(calc, 0, 6)[0]
    assert(uref.shape == (6, 6, 3))
    # NOTE: Comparing implementations of r gauge and v gauge is tricky, as they
    # tend to be numerically inequivalent.

    # compare to lrtddft implementation
    kss = KSSingles()
    atoms = calc.atoms
    atoms.calc = calc
    kss.calculate(calc.atoms, 1)
    lrref = []
    lrrefv = []
    for ex in kss:
        lrref.append(-1. * ex.mur * Bohr)
        lrrefv.append(-1. * ex.muv * Bohr)
    lrref = np.array(lrref)
    lrrefv = np.array(lrrefv)

    # Additional benefit: tests equivalence of r gauge implementations
    assert lrref[0, 2] == pytest.approx(uref[4, 0, 2])
    assert lrref[1, 1] == pytest.approx(uref[5, 0, 1])
    assert lrref[2, 1] == pytest.approx(uref[4, 1, 1])
    assert lrref[3, 2] == pytest.approx(uref[5, 1, 2])
    assert lrref[4, 2] == pytest.approx(uref[4, 2, 2])
    assert lrref[5, 1] == pytest.approx(uref[5, 2, 1])
    assert lrref[6, 0] == pytest.approx(uref[4, 3, 0])

    # some printout for manual inspection, if wanted
    print("         r-gauge   lrtddft(v)  raman(v)")
    f = "{} {:+.4f}    {:+.4f}    {:+.4f}"
    print(f.format('0->4 (z)', lrref[0, 2], lrrefv[0, 2], dip_vnm[2, 0, 4]))
    print(f.format('0->5 (y)', lrref[1, 1], lrrefv[1, 1], dip_vnm[1, 0, 5]))
    print(f.format('1->4 (y)', lrref[2, 1], lrrefv[2, 1], dip_vnm[1, 1, 4]))
    print(f.format('1->5 (z)', lrref[3, 2], lrrefv[3, 2], dip_vnm[2, 1, 5]))
    print(f.format('2->4 (z)', lrref[4, 2], lrrefv[4, 2], dip_vnm[2, 2, 4]))
    print(f.format('2->5 (y)', lrref[5, 1], lrrefv[5, 1], dip_vnm[1, 2, 5]))
    print(f.format('3->4 (x)', lrref[6, 0], lrrefv[6, 0], dip_vnm[0, 3, 4]))
