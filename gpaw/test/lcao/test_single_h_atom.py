import pytest
from ase import Atoms
from gpaw import GPAW, Mixer
# from gpaw.xc.noncolinear import NonColinearLDA, NonColinearLCAOEigensolver, \
#     NonColinearMixer


@pytest.mark.skip(reason='TODO')
def test_colinear():
    h = Atoms('H', magmoms=[1])
    h.center(vacuum=2)
    xc = 'LDA'
    c = GPAW(txt='c.txt',
             mode='lcao',
             basis='dz(dzp)',
             # setups='ncpp',
             h=0.25,
             xc=xc,
             # occupations=FermiDirac(0.01),
             nbands=1,
             mixer=Mixer(),
             # noncolinear=[(2,0,0)],
             )  # eigensolver=NonColinearLCAOEigensolver())
    h.calc = c
    h.get_potential_energy()
