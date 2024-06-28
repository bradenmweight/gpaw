.. _gbar submitting:

====================================
Submitting jobs on the DTU computers
====================================

Smaller calculations can be run in a Jupyter Notebook, but larger calculations
require running on multiple CPU cores for an extended time.  Such jobs should
be submitted with the MyQueue tool.  MyQueue_ is a unified frontend
for a number of different queuing systems available on HPC
installations.  It supports submitting individual jobs as well as
complete workflows.

.. _MyQueue: https://myqueue.readthedocs.io/en/latest/


Using MyQueue
=================

The command ``mq`` acts as a front-end to the queue system.
Usage::

  mq submit -R CORES:TIME script

Submit a GPAW Python script via the configured queueing system.

positional arguments:
  script:
    Python script

  argument:
    Command-line argument for Python script.

selected optional arguments:
  -h, --help            show help message and exit
  -n NAME, --name NAME  Name used for task.
  -R RESOURCES, --resources RESOURCES
                        Examples: "8:1h", 8 cores for 1 hour. Use "m" for minutes, "h" for hours and "d" for days. "16:1:30m": 16 cores,
                        1 process, half an hour.
  -z, --dry-run         Show what will happen without doing anything.
  -v, --verbose         More output.
  -q, --quiet           Less output.

.. code:: bash

    $ mq submit -R 8:4h script.py  # 8 cores, 4 hours
    $ mq list
    $ qstat hpc
    ...

The last command shows the user's jobs in the hpc queue, which is the
queue we use for the summer school.  ``mq list`` and ``qstat hpc`` give
some of the same information.


Choosing the number of processes
================================

GPAW parallelizes most efficiently over k-points, so it is a good idea to make
the number of processes a divisor of the number of *irreducible* k-points.  If
you have 12 irreducible k-points, the calculation parallelizes well on 2, 3,
4, 6 or 12 processes.

If you have very few irreducible k-points you may need to have more processes
than k-points; in these cases GPAW choose other parallelization strategies.
In this case, it is an advantage to make the number of processes a multiple of
the number of irreducible k-points.


Dry run: Let GPAW help you choosing
===================================

If you run your script with the command::

    $ gpaw python --dry-run=1 myscript.py

then your script will execute until the first GPAW calculation.  That
calculation will print information into the ``.txt`` file, and then stop.  In
the file, you can see the number of irreducible k-points and use it to select
your parallelization strategy.

Once you have decided how many processes you want, run another dry-run to
check how GPAW will parallelize::

    $ gpaw python --dry-run=PROCESSES myscript.py

where ``PROCESSES`` is the number of processes you want to use.  In this case,
gpaw will print how it will parallelize the calculation when running for real.
