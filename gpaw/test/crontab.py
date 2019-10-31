"""Test GPAW in a venv."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


cmds = """\
python3 -m venv venv
. venv/bin/activate
pip install wheel
pip install git+https://gitlab.com/ase/ase.git@master
pip install git+https://gitlab.com/gpaw/gpaw.git@master
gpaw test > test-1.out
gpaw -P 2 test > test-2.out
gpaw -P 4 test > test-4.out
gpaw -P 8 test > test-8.out"""


def run_tests():
    root = Path('/tmp/gpaw-tests')
    if root.is_dir():
        sys.exit('Locked')
    root.mkdir()
    os.chdir(root)
    cmds2 = ' && '.join(cmd for cmd in cmds.splitlines()
                        if not cmd.startswith('#'))
    p = subprocess.run(cmds2, shell=True)
    if p.returncode == 0:
        status = 'ok'
    else:
        print('FAILED!', file=sys.stdout)
        status = 'error'
    f = root.with_name('gpaw-test-' + status)
    if f.is_dir():
        shutil.rmtree(f)
    root.rename(f)
    return status


if __name__ == '__main__':
    run_tests()
