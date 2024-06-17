import pytest
from ase import Atoms
from ase.io import read

from gpaw import GPAW


@pytest.mark.ci
def test_no_cell():
    with pytest.raises(ValueError):
        H = Atoms('H', calculator=GPAW(mode='fd'))
        H.get_potential_energy()


@pytest.mark.parametrize('name', ['h2_pw', 'bcc_li_lcao'])
def test_read_txt(in_tmp_dir, gpw_files, name):
    gpw = gpw_files[name]
    e0 = GPAW(gpw).get_atoms().get_potential_energy()
    e = read(gpw.with_suffix('.txt')).get_potential_energy()
    assert e == pytest.approx(e0)
