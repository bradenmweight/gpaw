======
CentOS
======

Install these CentOS_ packages::

    $ yum install libxc-devel openblas-devel openmpi-devel fftw-devel
    $ yum install blacs-openmpi-devel scalapack-openmpi-devel

Add this to your ``~/.bashrc``::

    $ OPENMPI=/usr/lib64/openmpi
    $ export PATH=$OPENMPI/bin/:$PATH
    $ export LD_LIBRARY_PATH=$OPENMPI/lib:$LD_LIBRARY_PATH

Make sure you have the latest pip::

    $ python3 -m ensurepip --user
    $ python3 -m pip install pip --user

Then install ASE_, Numpy and SciPy::

    $ python3 -m pip install ase --user

And finally, GPAW with ScaLAPACK and FFTW::

    $ git clone git@gitlab.com:gpaw/gpaw.git
    $ cd gpaw
    $ cat > siteconfig.py
    fftw = True
    scalapack = True
    libraries = ['xc', 'fftw3', 'scalapack', 'mpiblacs']
    library_dirs = ['/usr/lib64/openmpi/lib/']
    ^D
    $ python3 -m pip install -v gpaw --user


.. _CentOS: https://www.centos.org/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
