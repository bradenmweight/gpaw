from myqueue.workflow import run
from pathlib import Path


def workflow():
    with run(script='H2O_ir.py'):
        run(script='H2O_ir_summary.py')
        with run(script='H2O_rraman_calc.py', cores=4, tmax='1h'):
            run(script='H2O_rraman_summary.py')
            run(script='H2O_rraman_spectrum.py')
            run(function=check)


def check():
    """Read result and make sure it's OK."""
    lines = Path('H2O_rraman_summary.txt').read_text().splitlines()
    for line in lines:
        if line.strip().startswith('8'):
            _, frequency, _, intensity = (float(x) for x in line.split())
            assert abs(frequency - 474) < 1
            assert abs(intensity - 57) < 1
            return
    raise ValueError
