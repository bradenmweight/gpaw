# Creates: Cu_bulk.csv


def workflow():
    from q2.job import Job
    return [
        Job('Cu_bulk.agts.py')]


if __name__ == '__main__':
    from ase.optimize.test.Cu_bulk import *
