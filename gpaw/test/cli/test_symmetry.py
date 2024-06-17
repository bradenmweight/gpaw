"""Test gpaw symmetry command."""
import pytest
from gpaw.cli.main import main

result = """\
symmetry:
  number of symmetries: 48
  number of symmetries with translation: 0

bz sampling:
  number of bz points: 512
  number of ibz points: 29
  monkhorst-pack size: [8, 8, 8]
  monkhorst-pack shift: [0.0625, 0.0625, 0.0625]

"""


@pytest.mark.serial
def test_symmetry(gpw_files, capsys):
    args = ['symmetry',
            str(gpw_files['bcc_li_pw']),
            '-k',
            '{density:3,gamma:1}']
    main(args)
    out = capsys.readouterr().out
    print(out)
    assert out == result
