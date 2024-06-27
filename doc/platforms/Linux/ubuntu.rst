=============
Ubuntu 18.04+
=============

.. highlight:: bash

Install these Ubuntu_ packages::

    $ sudo apt install python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev libfftw3-dev

Use this :download:`siteconfig.py` file:

.. literalinclude:: siteconfig.py

Put the file in your ``~/.gpaw/`` folder.  See also
:ref:`siteconfig`.

Then install GPAW (and dependencies: ASE_, Numpy, SciPy):

* Latest stable version from PyPI::

    $ pip install gpaw

* Development version:

  Clone the source code and install it::

    $ git clone https://gitlab.com/gpaw/gpaw.git
    $ pip install gpaw/


.. _Ubuntu: https://www.ubuntu.com/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
