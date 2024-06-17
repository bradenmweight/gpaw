.. _testing:

============
Testing GPAW
============

Testing of gpaw is done by a nightly test suite consisting of many
small and quick tests (with pytest) and by a weekly set of larger tests.


Test suite with pytest
======================

The test suite consists of a large number of small and quick tests
found in the :git:`gpaw/test/` directory.  The tests run nightly in serial
and in parallel modes.


Running tests in serial mode
----------------------------

Use pytest_ to run the tests::

    $ pytest --pyargs gpaw -v

To speed up the test suite, use pytest-xdist_ to use multiple processes
to run multiple tests at the same time
(note: each test is still run in serial mode)::

    $ pytest --pyargs gpaw -v -n <number-of-processes>

Please report errors to the ``gpaw-users`` mailing list so that we
can fix them (see :ref:`mail list`).

.. _pytest: http://doc.pytest.org/en/latest/contents.html
.. _pytest-xdist: https://github.com/pytest-dev/pytest-xdist


Running tests in parallel mode
------------------------------

In order to run the tests with MPI parallelization, do this::

    $ mpiexec -n <number-of-processes> pytest --pyargs gpaw -v

The tests should pass with 1, 2, 4, and 8 parallel tasks.

.. hint::

    If you observe issues (e.g. segmentation faults) when
    trying to run pytest, try this instead::

        $ mpiexec -n <n> gpaw python -m pytest --pyargs gpaw -v

    This should ensure that the correct environment is used.

Please report also parallel errors to the mailing list so that we
can fix them (see :ref:`mail list`).


Running a subset of tests
-------------------------

There are multiple options for running only a subset of test.

1. Use markers to run tests with that mark, for example CI tests::

    $ pytest --pyargs gpaw -v -m ci

2. Use module path to run tests in that path::

    $ pytest --pyargs gpaw.test.lcao -v

3. Use file/directory path to run tests in that path::

    $ pytest /root/of/gpaw/git/clone/gpaw/test/lcao


Special fixtures and marks
--------------------------

.. highlight:: python

Tests that should only run in serial can be marked like this::

    import pytest

    @pytest.mark.serial
    def test_something():
        ...

There are two special GPAW-fixtures:

.. autofunction:: gpaw.test.conftest.in_tmp_dir
.. autofunction:: gpaw.test.conftest.add_cwd_to_setup_paths
.. autofunction:: gpaw.test.conftest.gpw_files

Check the :git:`~gpaw/test/conftest.py` to see which gpw-files are available.
Use a ``_wfs`` post-fix to get a gpw-file that contains the wave functions.

.. autofunction:: gpaw.test.findpeak


Adding new tests
----------------

A test script should fulfill a number of requirements:

* It should be quick.  Preferably not more than a few milliseconds.
  If the test takes several minutes or more, consider making the
  test a :ref:`big test <big-test>`.

* It should not depend on other scripts.

* It should be possible to run it on 1, 2, 4, and 8 cores.

A test can produce standard output and files - it doesn't have to
clean up.  Just add the ``in_tmp_dir`` fixture as an argument::

    def test_something(in_tmp_dir):
        # make a mess ...

Here is a parametrized test that uses :func:`pytest.approx` for
comparing floating point numbers::

    import pytest

    @pytest.mark.parametrize('x', [1.0, 1.5, 2.0])
    def test_sqr(x):
        assert x**2 == pytest.approx(x * x)


.. _big-test:
.. _agts:

Big tests
=========

The directories in :git:`gpaw/test/big/` and :git:`doc/tutorialsexercises/`
contain longer and more
realistic tests that we run every weekend.  These are submitted to a
queuing system of a large computer.  The scripts in the :git:`doc` folder
are used both for testing GPAW and for generating up to date figures and
CSV-file for inclusion in the documentation web-pages.


Adding new tests
----------------

To add a new test, create a script somewhere in the file hierarchy ending with
``agts.py`` (e.g. ``submit.agts.py`` or just ``agts.py``). ``AGTS`` is short
for Advanced GPAW Test System (or Another Great Time Sink). This script
defines how a number of scripts should be submitted to Niflheim and how they
depend on each other. Consider an example where one script, ``calculate.py``,
calculates something and saves a ``.gpw`` file and another script,
``analyse.py``, analyses this output. Then the submit script should look
something like::

    def workflow():
        from myqueue.workflow import run
        with run(script='calculate.py', cores=8, tmax='25m'):
            run(script='analyse.py')  # 1 core and 10 minutes

As shown, this script has to contain the definition of the function
workflow_.  Start the workflow with ``mq workflow -p agts.py .``
(see https://myqueue.readthedocs.io/ for more details).

Scripts that generate figures or test files for inclusion in the
GPAW web-pages should start with a special ``# web-page:`` comment like this::

    # web-page: fig1.png, table1.csv
    ...
    # code that creates fig1.png and table1.csv
    ...

.. _workflow: https://myqueue.readthedocs.io/en/latest/
    workflows.html
