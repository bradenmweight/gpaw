# web-page: systems.db
from ase.optimize.test.test import all_optimizers
from ase.optimize.test.systems import create_database


def workflow():
    from myqueue.workflow import run

    with run(function=create_database):
        runs = [run(script='run_tests_emt.py')]

        for name in all_optimizers:
            if name in {'Berny', 'CellAwareBFGS'}:
                continue
            runs.append(run(script='run_tests.py',
                            args=[name], cores=8, tmax='1d'))

    run(script='analyze.py', deps=runs)
