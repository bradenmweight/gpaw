.. _devel:

===========
Development
===========

GPAW development can be done by anyone! Just take a look at the
:ref:`todo` list and find something that suits your talents!

The primary source of information is still the :ref:`manual` and
:ref:`documentation`, but as a developer you might need additional
information which can be found here. For example the :ref:`code_overview`.

As a developer, you should subscribe to all GPAW related :ref:`mailing_lists`.

Now you are ready to to perfom a :ref:`developer_installation` and
start development!

.. toctree::
   :maxdepth: 1

   developer_installation

.. note --- below toctrees are defined in separate files to make sure that the line spacing doesn't get very large (which is of course a bad hack)

Development topics
==================

When committing significant changes to the code, remember to add a
note in the :ref:`releasenotes` at the top (current svn) - the version
to become the next release.

.. toctree::
   :maxdepth: 1

   toc-general

* The latest report_ from PyLint_ about the GPAW coding standard.

.. spacer

* Details about supported :ref:`platforms_and_architectures`.

.. _report: http://dcwww.camd.dtu.dk/~s052580/pylint/gpaw/
.. _PyLint: http://www.logilab.org/857


.. _code_overview:

Code Overview
=============

Keep this picture under your pillow:

.. _the_big_picture:

.. image:: bigpicture.png
   :target: ../bigpicture.pdf

The developer guide provides an overview of the PAW quantities and how
the corresponding objects are defined in the code:

.. toctree::
   :maxdepth: 2

   overview
   developersguide
   paw
   wavefunctions
   setups
   density_and_hamiltonian
   communicators
   others


The GPAW logo
=============

The GPAW-logo is available in the odg_ and svg_ file formats:
gpaw-logo.odg_, gpaw-logo.svg_

.. _odg: http://www.openoffice.org/product/draw.html
.. _svg: http://en.wikipedia.org/wiki/Scalable_Vector_Graphics
.. _gpaw-logo.odg: ../_static/gpaw-logo.odg
.. _gpaw-logo.svg: ../_static/gpaw-logo.svg


Statistics
==========

The image below shows the development in the volume of the code as per
December 2 2009.

.. image:: ../_static/stat.png

*Documentation* refers solely the contents of this homepage. Inline
documentation is included in the other line counts.

.. ::

   Commented text!

   The gpaw development project currently count about 20 active
   developers in 7 universities, and has 93 subscribers to the mailing
   list.


