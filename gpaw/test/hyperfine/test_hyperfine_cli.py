"""Test hyperfine parameter CLI."""
import pytest
import numpy as np

from gpaw.hyperfine import main


@pytest.mark.parametrize('option', ['', '-x', '-d'])
@pytest.mark.serial
def test_hyperfine_cli(gpw_files, capsys, option):
    file = gpw_files['o2_pw']
    main(f'{file} -g O:-0.757 -u MHz {option}'.split())
    out = capsys.readouterr().out
    print(out)
    lines = [line.split() for line in out.splitlines()]
    if option == '-x':
        ref = [0.689, 33.74, -246.58, -246.58, 0.0, 0.0, 0.0]
    elif option == '-d':
        ref = [0.689, -78.01, -78.01, 202.31]
    else:
        ref = [0.689, 202.31, -78.01, -78.01, 0.0, 0.0, 0.0]
    for line in lines[3:5]:
        errors = (np.array([float(x) for x in line[2:]]) - ref)
        assert abs(errors).max() < 0.2
    assert float(lines[7][3]) == 2.0
    assert float(lines[10][1]) == -0.757
