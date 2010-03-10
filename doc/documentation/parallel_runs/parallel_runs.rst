.. _parallel_runs:

=============
Parallel runs
=============

.. toctree::
   :maxdepth: 1

.. _parallel_running_jobs:

Running jobs in parallel
========================

Parallel calculations are done with MPI and a special
:program:`gpaw-python` python-interpreter.

The parallelization can be done over the **k**-points, bands, spin in
spin-polarized calculations, and using real-space domain
decomposition.  The code will try to make a sensible domain
decomposition that match both the number of processors and the size of
the unit cell.  This choice can be overruled, see
:ref:`manual_parallelization_types`.

Before starting a parallel calculation, it might be useful to check how the parallelization corresponding to given number
of processors would be done with ``--dry-run`` command line option::

  python script.py --dry-run=8

In order to start parallel calculation, you need to know the
command for starting parallel processes. This command might contain
also the number of processors to use and a file containing the names
of the computing nodes.  Some
examples::

  mpirun -np 4 gpaw-python script.py
  poe "gpaw-python script.py" -procs 8


Simple submit tool
==================

Instead writing a file with the line "mpirun ... gpaw-python script.py" and then submitting it to a queueing system, it is simpler to automate this::

  #!/usr/bin/env python
  from sys import argv
  import os
  options = ' '.join(argv[1:-1])
  job = argv[-1]
  dir = os.getcwd()
  f = open('script.sh', 'w')
  f.write("""\
  NP=`wc -l < $PBS_NODEFILE`
  cd %s
  mpirun -np $NP -machinefile $PBS_NODEFILE gpaw-python %s
  """ % (dir, job))
  f.close()
  os.system('qsub ' + options + ' script.sh')

Now you can do::

  $ qsub.py -l nodes=20 -m abe job.py

You will have to modify the script so that it works with your queueing
system.

.. _submit_tool_on_niflheim:

Submit tool on Niflheim
-----------------------

At CAMd, we use this submit tool: :svn:`~doc/documentation/parallel_runs/gpaw-qsub`.

Example::

  $ gpaw-qsub -q medium -l nodes=8 -m abe fcc.py --domain-decomposition=1,2,2

.. tip::
  CAMd users must always remember to source the openmpi environment settings before recompiling the code. See :ref:`Niflheim`.

Alternative submit tool
=======================

Alternatively, the script gpaw-runscript can be used, try::

  $ gpaw-runscript -h

to get the architectures implemented and the available options. As an example, use::

  $ gpaw-runscript script.py 32

to write a job sumission script running script.py on 32 cpus. The tool tries to guess the architecture/host automatically.


Writing to files
================

Be careful when writing to files in a parallel run.  Instead of ``f = open('data', 'w')``, use:

>>> from ase.parallel import paropen
>>> f = paropen('data', 'w')

Using ``paropen``, you get a real file object on the master node, and dummy objects on the slaves.  It is equivalent to this:

>>> from ase.parallel import rank
>>> if rank == 0:
...     f = open('data', 'w')
... else:
...     f = open('/dev/null', 'w')

If you *really* want all nodes to write something to files, you should make sure that the files have different names:

>>> from ase.parallel import rank
>>> f = open('data.%d' % rank, 'w')

Writing text output
===================

Text output written by the ``print`` statement is written by all nodes.
To avoid this use:

>>> from ase.parallel import parprint
>>> print 'This is written by all nodes'
>>> parprint('This is written by the master only')

which is equivalent to

>>> from ase.parallel import rank
>>> print 'This is written by all nodes'
>>> if rank == 0:
...     print 'This is written by the master only'

Note that parprint has to be used as the print statement of python 3.

Running different calculations in parallel
==========================================

A GPAW calculator object will per default distribute its work on all
available processes. If you want to use several different calculators
at the same time, however, you can specify a set of processes to be
used by each calculator. The processes are supplied to the constructor,
either by specifying an :ref:`MPI Communicator object <communicators>`,
or simply a list of ranks. Thus, you may write::

  from gpaw import GPAW
  import gpaw.mpi as mpi
  import numpy as np

  # Create a calculator using ranks 0, 3 and 4 from the mpi world communicator
  comm = mpi.world.new_communicator(np.array([0, 3, 4]))
  calc = GPAW(communicator=comm)

Be sure to specify different output files to each calculator,
otherwise their outputs will be mixed up.

Here is an example which calculates the atomization energy of a
nitrogen molecule using two processes:

.. literalinclude:: parallel_atomization.py

ScaLapack
=========

.. toctree::
   :maxdepth: 1

   ScaLapack/ScaLapack

.. _manual_parallelization_types:

Parallization modes
===================

.. _manual_parsize:

Domain decomposition
--------------------

The choice for the domain decomposition can be forced using the 
keyword ``parsize``. It can be set like ``parsize=(nx,ny,nz)`` to
force the decomposition into nx, ny, and nz boxes in x,y, and z direction
respectively. The use of domain decomposition only can be forced by
``parsize='domain only'``.

There is also a command line argument
``--domain-decomposition`` that allows you to control the domain
decomposition (see example at :ref:`submit_tool_on_niflheim`).

Band parallelization
--------------------

.. toctree::
   :maxdepth: 1

   band_parallelization/band_parallelization
