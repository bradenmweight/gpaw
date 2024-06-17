.. _installation:

============
Installation
============

.. toctree::
    :hidden:

    troubleshooting
    platforms/platforms

GPAW relies on the Python library *atomic simulation environment* (ASE_),
so you need to :ref:`install ASE <ase:download_and_install>` first.  GPAW
itself is written mostly in the Python programming language, but there
are also some C-code used for:

* performance critical parts
* allowing Python to talk to external numerical libraries (BLAS_,
  LibXC_, MPI_ and ScaLAPACK_)

So, in order to make GPAW work, you need to compile some C-code.
You will need to build a dynamically linked library
(``_gpaw.so``) that the standard Python interpreter can load.

There are several ways to install GPAW:

* Directly from PyPI_: ``python3 -m pip install gpaw``.  More details
  here: :ref:`installation using pip`.

* Alternatively, you can :ref:`download <download>` the source code
  and then install with a
  ``python3 -m pip install .``.

* There may be a package for your Linux distribution that you can use
  (named ``gpaw``).

* If you are a developer that need to change the code you should look at this
  description: :ref:`development workflow`.

.. seealso::

    * :ref:`siteconfig`
    * Using :ref:`homebrew` on MacOSX
    * Using :ref:`anaconda`
    * This `docker image`_
    * Tips and tricks for installation on many :ref:`platforms and
      architectures`
    * :ref:`troubleshooting`
    * In case of trouble: :ref:`Our mail list and #gpaw chat channel <contact>`


Requirements
============

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

See :ref:`releasenotes` for version requirements.

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _LibXC: http://www.tddft.org/programs/libxc/
.. _MPI: http://www.mpi-forum.org/
.. _BLAS: http://www.netlib.org/blas/
.. _BLACS: http://www.netlib.org/blacs/
.. _ScaLAPACK: http://www.netlib.org/scalapack/
.. _PyPI: https://pypi.org/project/gpaw
.. _PIP: https://pip.pypa.io/en/stable/
.. _ASE: https://wiki.fysik.dtu.dk/ase
.. _FFTW: http://www.fftw.org/
.. _docker image: https://hub.docker.com/r/marcindulak/gpaw-openmpi


.. _installation using pip:

Installation using ``pip``
==========================

.. highlight:: bash

The simplest way to install GPAW is using pip_ and the GPAW package from
the Python package index (PyPI_)::

    $ python3 -m pip install gpaw

This will compile and install GPAW (both ``_gpaw.so`` and all the Python
files) in your ``~/.local/lib/pythonX.Y/site-packages`` folder where
Python can automatically find it.  Pip will also place
the command line tool :command:`gpaw` in the ``~/.local/bin`` folder, so
make sure you have that in your ``$PATH`` environment variable.

Check that you have installed everything in the correct places::

    $ gpaw info


Install PAW datasets
====================

Install the datasets into the folder ``<dir>`` using this command::

    $ gpaw install-data <dir>

See :ref:`installation of paw datasets` for more details.

Now you should be ready to use GPAW, but before you start, please run the
tests as described below.


.. index:: test
.. _run a tests:

Run a simple test calculation
=============================

Make sure that everything works by running a simple test::

    $ gpaw test

If this worked then try also a parallel calculation::

    $ gpaw -P 4 test

Please report errors to the ``gpaw-users`` mailing list so that we
can fix them (see :ref:`mail list`).

.. seealso::

    :ref:`testing` for how to run the complete test suite.


.. _download:

Getting the source code
=======================

You can get the source from a tar-file or from Git:

:Tar-file:

    You can get the source as a tar-file for the
    latest stable release (gpaw-24.6.0.tar.gz_) or the latest
    development snapshot (`<snapshot.tar.gz>`_).

    Unpack and make a soft link::

        $ tar -xf gpaw-24.6.0.tar.gz
        $ ln -s gpaw-24.6.0 gpaw

    Here is a `list of tarballs <https://pypi.org/simple/gpaw/>`__.

:Git clone:

    Alternatively, you can get the source for the latest stable release from
    https://gitlab.com/gpaw/gpaw like this::

        $ git clone -b 24.6.0 https://gitlab.com/gpaw/gpaw.git

    or if you want the development version::

        $ git clone https://gitlab.com/gpaw/gpaw.git

.. note::

    We also have Git tags for older stable versions of GPAW.
    See the :ref:`releasenotes` for which tags are available.  Also the
    dates of older releases can be found there.

.. _gpaw-24.6.0.tar.gz:
    https://pypi.org/packages/source/g/gpaw/gpaw-24.6.0.tar.gz


.. _siteconfig:

Customizing installation
========================

The install script may need a little help finding you libraries
(BLAS, FFTW, ScaLapack, libxc, libvdwxc, ...).
This can be done by adding a ``siteconfig.py`` file in one of these three
places:

1) the file that ``$GPAW_CONFIG`` points at
2) ``<git-root>/siteconfig.py``
3) ``~/.gpaw/siteconfig.py``

The first one found will be used. As an example, ``siteconfig.py``
might contain the following lines::

  libraries = ['myblas']
  library_dirs = ['path_to_myblas']

Now, GPAW would be built with "``-Lpath_to_myblas -lmyblas``"
linker flags. Look at the file :git:`siteconfig_example.py`
for more possible options.  :ref:`platforms and architectures`
provides examples of ``siteconfig.py`` files for different platforms.


.. _parallel installation:

Parallel installation
=====================

By default, setup looks if :program:`mpicc` is available, and if setup
finds one, a parallel version is build. If the setup does not find
:program:`mpicc`, a user can specify one in the ``siteconfig.py`` file.

Additionally a user may want to enable ScaLAPACK, setting in
``siteconfig.py``::

    scalapack = True

and, in this case, provide BLACS/ScaLAPACK ``libraries`` and ``library_dirs``
as described in :ref:`siteconfig`.

Part of the code can be parallelized also with OpenMP. OpenMP
threading can be enabled by setting in ``siteconfig.py``::

    extra_compile_args += ['-fopenmp']
    extra_link_args += ['-fopenmp']

Note that a multithreaded BLAS library should be used together with
OpenMP parallelization. Also, MPI library needs to support
the MPI_THREAD_MULTIPLE thread support level.

Instructions for running parallel calculations can be found in the
:ref:`user manual <parallel_runs>`.


FFTW
====

In order to use FFTW_ instead of :mod:`numpy.fft`, add something like
this to your ``siteconfig.py``::

    fftw = True
    libraries += ['fftw3']



.. _libxc installation:

Libxc Installation
==================

If you OS does not have a LibXC package you can use then you can download
and install LibXC as described `here <http://www.tddft.org/programs/libxc/>`_.
A few extra tips:

* Libxc installation requires both a C compiler and a Fortran compiler.

* We've tried Intel and gnu compilers and haven't noticed much of a
  performance difference.  Use whatever is easiest.

* Libxc shared libraries can be built with the "--enable-shared" option
  to configure.  This might be slightly preferred because it reduces
  memory footprints for executables.

* Typically when building GPAW one has to modify ``siteconfig.py`` in a manner
  similar to the following::

    library_dirs += ['/my/path/to/libxc/5.2.0/install/lib']
    include_dirs += ['/my/path/to/libxc/5.2.0/install/include']

  or if you don't want to modify your ``siteconfig.py``, you can add these
  lines to your ``.bashrc``::

    export C_INCLUDE_PATH=/my/path/to/libxc/5.2.0/install/include
    export LIBRARY_PATH=/my/path/to/libxc/5.2.0/install/lib
    export LD_LIBRARY_PATH=/my/path/to/libxc/5.2.0/install/lib

Example::

    wget http://www.tddft.org/programs/octopus/down.php?file=libxc/5.2.0/libxc-5.2.0.tar.gz -O libxc-5.2.0.tar.gz
    tar -xf libxc-5.2.0.tar.gz
    cd libxc-5.2.0
    ./configure --enable-shared --disable-fortran --prefix=$HOME/libxc-5.2.0
    make
    make install

    # add these to your .bashrc:
    XC=~/libxc-5.2.0
    export C_INCLUDE_PATH=$XC/include
    export LIBRARY_PATH=$XC/lib
    export LD_LIBRARY_PATH=$XC/lib
