import pytest
from gpaw.test.fuzz import main


@pytest.mark.ci
@pytest.mark.parametrize('mode', ['pw', 'lcao', 'fd'])
@pytest.mark.parametrize('pbc', [0, 1])
def test_fuzz(tmp_path, mode, pbc):
    error = main(f'h2 -m {mode} -c new,old -v 3.0 -p {pbc}')
    assert error == 0, mode
