.. _libvdwxc-doc:

libvdwxc
========

`libvdwxc <https://gitlab.com/libvdwxc/libvdwxc>`_
is a library which provides fast and scalable
implementations of non-local van der Waals density functionals in the
vdW-DF family.
To use libvdwxc, you need to install it
and compile GPAW with it.  libvdwxc can be used with other semilocal
functionals like optPBE, optB88, and BEEF-vdW.

`Install <https://libvdwxc.readthedocs.io>`_ libvdwxc,
making sure that its dependencies FFTW3 and
FFTW3-MPI are available on the system.  For truly large systems, you
may install PFFT to achieve better scalability.  For realistically-sized systems, FFTW3-MPI is efficient
and might be a bit faster than PFFT.

Run a calculation by specifying backend, like {'name':'BEEF-vdW', 'backend':'libvdwxc'}, as in this example:

.. literalinclude:: libvdwxc-example.py

libvdwxc will automatically parallelize with as many cores as are
available for domain decomposition.  If you parallelize over *k*-points
or bands, and *especially* if you use planewave mode, be sure to pass
the parallelization keyword ``augment_grids=True`` to make use of *all*
cores including those for *k*-point and band parallelization (see :ref:`parallel_runs`).

Note that libvdwxc 0.4 has no stress term implementation.
