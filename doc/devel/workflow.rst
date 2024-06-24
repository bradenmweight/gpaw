.. _development workflow:

====================
Development workflow
====================

.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _NumPy: https://docs.scipy.org/doc/numpy/reference/
.. _SciPy: https://docs.scipy.org/doc/scipy/reference/
.. _venv: https://docs.python.org/3/library/venv.html#module-venv
.. _pip: https://pip.pypa.io/
.. _git: https://git-scm.com/
.. _GitLab issues: https://gitlab.com/gpaw/gpaw/issues
.. _pytest: https://docs.pytest.org/en/6.2.x/

.. contents::

.. seealso::

   * :ref:`writing documentation`
   * :ref:`testing`


Setting up your development environment
=======================================

Make a `virtual environment <venv>`_::

 $ mkdir devel
 $ cd devel
 $ unset PYTHONPATH
 $ python3 -m venv venv
 $ source venv/bin/activate  # venv/bin/ is now first in $PATH
 $ pip install --upgrade pip

Install master branch of ASE_ in *editable* mode::

 $ git clone git@gitlab.com:ase/ase
 $ pip install --editable ase/

Same thing for GPAW::

 $ git clone git@gitlab.com:gpaw/gpaw
 $ echo "noblas = True; nolibxc = True" > gpaw/siteconfig.py
 $ pip install -e gpaw

.. note::

    Here we used a simple ``siteconfig.py`` that *should* always work:

    * ``noblas = True``: Use the BLAS library built into  NumPy_
      (usually OpenBLAS).
    * ``nolibxc = True``: Use GPAW's own XC-functionals
      (only LDA, PBE, revPBE, RPBE and PW91).

    See :ref:`siteconfig` for details.

Download PAW datasets::

 $ gpaw install-data --register ~/PAWDATA


Run the tests
=============

The test-suite can be found in :git:`gpaw/test/`.  Run it like this::

 $ pip install pytest-xdist
 $ cd gpaw
 $ pytest -n4

And with MPI (2, 4 and 8 cores)::

 $ mpiexec -n 2 pytest

.. warning::

   This will take forever!  It's a good idea to learn and master pytest_'s
   command-line options for selecting the subset of all the tests that are
   relevant.


Creating a merge request
========================

Request to become a member of the ``gpaw`` project on GitLab
`here <https://gitlab.com/gpaw/gpaw/>`__.  This will
allow you to push branches to the central repository (see below).

Create a branch for your changes::

 $ cd gpaw
 $ git switch -c fix-something

.. note::

   ``git switch -c fix-something`` is the same as any of these:

   * ``git branch fix-something && git switch fix-something``
   * ``git branch fix-something && git checkout fix-something``
   * ``git checkout -b fix-something``

   :xkcd:`More git-tricks <1597>`.

Make some changes and commit::

 $ git add file1 file2 ...
 $ git commit -m "Short summary of changes"

Push your branch to GitLab::

 $ git push --set-upstream origin fix-something

and click the link to create a merge-request (MR).  Mark the MR as DRAFT to
signal that it is work-in-progress and remove the DRAFT-marker once the MR
is ready for code review.

Every time you push your local repository changes upstream to the remote
repository, you will trigger a continuous integration (CI) runner on the
GitLab servers.  The script that runs in CI is :git:`.gitlab-ci.yml`.
Here is a short summary of what happens in CI:

* install the code
* ``pytest -m ci``: small selection of fast tests
* ``mypy -p gpaw``: `Static code analysis`_ (type hints)
* ``flake8``: pyflakes + pycodestyle (pep8) = flake8_

If CI fails, you will have to fix things and push your changes.

It's a good idea to also run the CI-checks locally::

 $ pip install flake8 mypy
 $ flake8 ...
 $ mypy ...
 $ pytest ...
 $ # fix things
 $ git add ...
 $ git commit ...
 $ git push  # Git now knows your upstream

.. tip::

   You can use ``git push -i ci.skip`` if you want to skip CI.


.. _Static code analysis: https://mypy.readthedocs.io/en/stable/
.. _flake8: https://flake8.pycqa.org/en/latest/


How to write a good MR
======================

A good MR

* is short
* does one thing
* is not too old

For MRs with code changes:

* make sure there is a test that covers the new/fixed code
* make sure all variable and functions have descriptive names.
* remember docstrings - if needed
  (no need for an ``add_numbers()`` function to have an
  ``"""Add numbers."""`` docstring).

For MRs with documentation changes,
build the HTML-pages and make sure everything looks OK::

 $ pip install sphinx-rtd-theme
 $ cd gpaw/doc
 $ make
 $ make browse
