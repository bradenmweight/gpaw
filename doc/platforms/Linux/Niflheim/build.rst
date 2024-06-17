.. _build on niflheim:

==========================================
Building GPAW in a Python venv on Niflheim
==========================================

This document explains how to compile a developer version of GPAW on
Niflheim.  If you just want to run the pre-installed version, please
read the guide :ref:`Using a pre-installed GPAW on Niflheim <load on niflheim>`.

.. seealso::

    * :mod:`Creation of Python virtual environments <venv>`.
    * Information about the Niflheim cluster can be found at
      `<https://niflheim-users.readthedocs.io/>`_ and
      `<https://niflheim-system.readthedocs.io/>`_.
    * `MyQueue <https://myqueue.readthedocs.io/>`__.

.. contents::

.. highlight:: bash


Creating the venv
=================

Download the :download:`gpaw_venv.py` script (use ``wget <url>``)
and run it like this::

    $ python3 gpaw_venv.py <venv-name>
    ...

.. tip::

    You will need Python 3.8 or later.  You can install that with::

        $ module load Python

Type ``python3 gpaw_venv.py --help`` for help.  After a few minutes, you will
have a ``<venv-name>`` folder with a GPAW installation inside.

In the following, we will assume that your venv folder is ``~/venv1/``.

The ``gpaw_venv.py`` script does the following:

* load relevant modules from the foss toolchain
* create the venv
* clone and install ASE and GPAW from gitlab
* install some other Python packages from PyPI: sklearn, graphviz,
  matplotlib, pytest-xdist, myqueue, ase-ext, spglib
* enable tab-completion for command-line tools:
  `ase <https://wiki.fysik.dtu.dk/ase/cmdline.html>`__,
  `gpaw <https://wiki.fysik.dtu.dk/gpaw/documentation/cmdline.html>`__,
  `mq <https://myqueue.readthedocs.io/en/latest/cli.html>`__


Using the venv
==============

The venv needs to be activated like this::

    $ source venv1/bin/activate

and you can deactivate it when you no longer need to use it::

    $ deactivate

You will want the activation to happen automatically for the jobs you
submit to Niflheim.  Here are three ways to do it (pick one, and only one):

1) If you always want to use one venv then just put the activation
   command in your ``~/.bashrc``.

2) If you only want jobs running inside a certain folder to use the venv,
   then add this to your ``~/.bashrc``::

       if [[ $SLURM_SUBMIT_DIR/ = $HOME/project-1* ]]; then
           source ~/venv1/bin/activate
       fi

   Now, SLURM-jobs submitted inside your ``~/project-1/``
   folder will use the venv.

3) Use MyQueue.  Make sure you have MyQueue version 22.7.0 or later
   (``mq --version``).  The venv will automatically be activated if it was
   activated at submit time.

   If you haven't configured MyQueue then you can do that with this command::

       $ mq config slurm | grep -v sm3090 > ~/.myqueue/config.py

   (skips the *sm3090* GPU-enabled nodes).



Adding additional packages
==========================

In order to add more Python packages to your venv, you need to activate it
and then you can ``pip install`` packages.  Here is how
to install ASR_::

    $ git clone https://gitlab.com/asr-dev/asr.git
    $ cd asr
    $ git checkout old-master
    $ pip install .

.. warning::

    Pip may need co compile some code.
    It is therefore safest to use the ``thul`` login node to pip install
    software as it is the oldest CPU architcture and the other login nodes
    will understand its code.

.. _ASR: https://asr.readthedocs.io/en/latest/


Full script
===========

.. literalinclude:: gpaw_venv.py
