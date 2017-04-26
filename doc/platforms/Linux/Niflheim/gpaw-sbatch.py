#!/usr/bin/env python

import sys
import subprocess


def main():
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg.endswith('.py'):
            break
    else:
        print('Submit a GPAW Python script via sbatch.')
        print('Usage: gpaw-sbatch [sbatch options] script.py '
              '[script arguments]')
        return

    mpi_cmd = 'mpiexec "$GPAW_MPI_OPTIONS" '
    mpi_cmd += 'gpaw-python '
    script = '#!/bin/bash -l\n'
    for line in open(arg):
        if line.startswith('#SBATCH'):
            script += line
    script += ('OMP_NUM_THREADS=1 ' +
               mpi_cmd +
               ' '.join(args[i:]) + '\n')
    cmd = ['sbatch'] + args[:i]
    subprocess.Popen(cmd, stdin=subprocess.PIPE).communicate(script.encode())


if __name__ == '__main__':
    main()
