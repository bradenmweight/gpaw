.. _devel:

===========
Development
===========

To get started, take a look a our :ref:`development workflow`.

GPAW development can be done by anyone! Just take a look at the
`issue tracker`_ and find something that suits your talents.

The primary source of information is still the :ref:`basic` and
:ref:`documentation`, but as a developer you might need additional
information which can be found here. For example the :ref:`code_overview`.

As a developer, you should subscribe to the GPAW :ref:`mail list`.
We would also like to encourage you to join our ``#gpaw`` channel
for :ref:`chat`.

.. _issue tracker: https://gitlab.com/gpaw/gpaw/issues/

.. toctree::
   :maxdepth: 2

   workflow

.. note --- below toctrees are defined in separate files to make sure that the line spacing doesn't get very large (which is of course a bad hack)


Development topics
==================

When committing significant changes to the code, remember to add a
note in the :ref:`releasenotes` at the top (development version) - the
version to become the next release.

.. toctree::
   :maxdepth: 1

   testing
   codingstandard
   c_extension
   writing_documentation
   formulas
   debugging
   turn_off_things
   profiling
   ase_optimize/ase_optimize
   bugs
   newrelease
   technology

* Details about supported :ref:`platforms and architectures`.

.. _PyLint: https://www.logilab.org/857


.. _the_big_picture:
.. _code_overview:

Code Overview
=============

The developer guide provides an overview of the PAW quantities and how
the corresponding objects are defined in the code:

.. toctree::
   :maxdepth: 2

   overview
   developersguide
   proposals/proposals
   paw
   symmetry
   wavefunctions
   setups
   density_and_hamiltonian
   communicators
   others


The GPAW logo
=============

The GPAW-logo is available as an SVG-file: :download:`gpaw-logo.svg`.

.. image:: gpaw-logo.svg


Statistics
==========

The image below shows the development in the volume of the code as per
April 5 2016.

.. image:: lines.png

*Documentation* refers solely the contents of this homepage. Inline
documentation is included in the other line counts.


Contributing to GPAW
====================

Getting commit access to the GPAW code works the same way as for
the :ref:`ASE project <ase:contribute>`.
