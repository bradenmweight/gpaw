.. image:: https://badge.fury.io/py/gpaw.svg
    :target: https://pypi.org/project/gpaw/

GPAW
====

GPAW is a density-functional theory (DFT) Python_ code based on the
projector-augmented wave (PAW) method and the atomic simulation environment
(ASE_). It uses plane-waves, atom-centered basis-functions or real-space
uniform grids combined with multigrid methods.

Webpage: https://wiki.fysik.dtu.dk/gpaw


Requirements
------------

* Python_ 3.8 or later
* ASE_
* NumPy_
* SciPy_
* LibXC_
* A C compiler
* A BLAS_ library

Optional (highly recommended for increased performance):

* A MPI_ library (required for parallel calculations)
* BLACS_ and ScaLAPACK_ libraries
* FFTW_

See `Release notes <releasenotes_>`_ for version requirements.


Installation
------------

Create a virtual environment, activate it, install::

    $ python3 -m venv venv
    $ source venv/bin/activare
    $ python3 -m pip install gpaw

For more details, please see:

    https://wiki.fysik.dtu.dk/gpaw/install.html


Test your installation
----------------------

You can do a test calculation with::

    $ gpaw test


Contact
-------

* Mailing list: gpaw-users_
* Chat: #gpaw on Matrix_.
* Bug reports and development: gitlab-issues_

Please send us bug-reports, patches, code, ideas and questions.


Example
-------

Geometry optimization of hydrogen molecule:

>>> from ase import Atoms
>>> from ase.optimize import BFGS
>>> from ase.io import write
>>> from gpaw import GPAW, PW
>>> h2 = Atoms('H2',
...            positions=[[0, 0, 0],
...                       [0, 0, 0.7]])
>>> h2.center(vacuum=2.5)
>>> h2.calc = GPAW(xc='PBE',
...                mode=PW(300),
...                txt='h2.txt')
>>> opt = BFGS(h2, trajectory='h2.traj')
>>> opt.run(fmax=0.02)
BFGS:   0  09:08:09       -6.566505       2.2970
BFGS:   1  09:08:11       -6.629859       0.1871
BFGS:   2  09:08:12       -6.630410       0.0350
BFGS:   3  09:08:13       -6.630429       0.0003
>>> write('H2.xyz', h2)
>>> h2.get_potential_energy()  # ASE's units are eV and Å
-6.6304292169392784


Getting started
---------------

Once you have familiarized yourself with ASE_ and NumPy_, you should take a
look at the GPAW exercises_ and tutorials_.


.. _Python: https://www.python.org/
.. _ASE: https://wiki.fysik.dtu.dk/ase
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _LibXC: https://libxc.gitlab.io/
.. _MPI: http://www.mpi-forum.org/
.. _BLAS: http://www.netlib.org/blas/
.. _FFTW: http://www.fftw.org/
.. _BLACS: http://www.netlib.org/blacs/
.. _ScaLAPACK: http://www.netlib.org/scalapack/
.. _gpaw-users: https://listserv.fysik.dtu.dk/mailman/listinfo/gpaw-users
.. _Matrix: https://app.element.io/#/room/#gpaw:matrix.org
.. _gitlab-issues: https://gitlab.com/gpaw/gpaw/issues
.. _exercises: https://wiki.fysik.dtu.dk/gpaw/exercises/exercises.html
.. _tutorials: https://wiki.fysik.dtu.dk/gpaw/tutorials/tutorials.html
.. _releasenotes: https://wiki.fysik.dtu.dk/gpaw/releasenotes.html
