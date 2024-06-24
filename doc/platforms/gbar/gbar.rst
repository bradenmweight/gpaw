========
DTU GBar
========

This document shows how to build a *venv* for GPAW+ASE.

.. seealso::

    * :mod:`Creation of Python virtual environments <venv>`.
    * Information about the `GBar <https://www.gbar.dtu.dk/>`_.
    * Information about the `DCC software stack
      <https://www.hpc.dtu.dk/?page_id=3246>`_.
    * MyQueue_.


.. _MyQueue: https://myqueue.readthedocs.io/
.. highlight:: bash


Creating the venv
=================

Download the ``gpaw-venv.sh`` script
using this link: :download:`gpaw-venv.sh` or these commands::

    $ gpaw=https://gitlab.com/gpaw/gpaw
    $ wget $gpaw/-/raw/master/doc/platforms/gbar/gpaw-venv.sh

and run it like this::

    $ sh gpaw-venv.sh <venv-name>
    ...

After a few minutes, you will have a ``<venv-name>`` folder with
a GPAW installation inside (plus some other stuff).

.. note::

    The GPAW installation will only work on the
    ``XeonE5_2650v4`` and ``XeonE5_2660v3`` architectures.

In the following, we will assume that your venv folder is ``~/venv/``.


Using the venv
==============

The venv needs to be activated like this::

    $ source ~/venv/bin/activate

and you can deactivate it when you no longer need to use it::

    $ deactivate


Submitting jobs
===============

Using bsub
----------

See `here <https://www.gbar.dtu.dk/>`_.


Using MyQueue
-------------

First, configure MyQueue_::

    $ mq config --in-place -Q hpc lsf

.. warning::

   Your ``~/.myqueue/config.py`` file should look like this:

   .. literalinclude:: config.py

   Edit this file so that it only contains the
   ``XeonE5_2650v4`` and ``XeonE5_2660v3`` architectures.

Then you can submit jobs with::

    $ mq submit script.py -R8:4h  # 8 cores, 5 hours
    $ mq ls
    $ mq --help
