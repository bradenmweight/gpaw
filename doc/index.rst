===============================================================
GPAW: DFT and beyond within the projector-augmented wave method
===============================================================

GPAW is a density-functional theory (DFT) Python_ code based on the
projector-augmented wave (:ref:`PAW <introduction_to_paw>`) method and the
atomic simulation environment (ASE_).  The wave functions can be described
with:

* Plane-waves (:ref:`pw <manual_mode>`)
* Real-space uniform grids, multigrid methods and the finite-difference approximation
  (:ref:`fd <manual_stencils>`)
* Atom-centered basis-functions (:ref:`lcao <lcao>`)

>>> # H2-molecule example:
>>> import numpy as np
>>> from ase import Atoms
>>> from gpaw import GPAW, PW
>>> h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.74)])
>>> h2.center(vacuum=2.5)
>>> h2.cell
Cell([5.0, 5.0, 5.74])
>>> h2.positions
array([[2.5 , 2.5 , 2.5 ],
       [2.5 , 2.5 , 3.24]])
>>> h2.calc = GPAW(xc='PBE',
...                mode=PW(300),
...                txt='h2.txt')
>>> energy = h2.get_potential_energy()
>>> print(f'Energy: {energy:.3f} eV')
Energy: -6.631 eV
>>> forces = h2.get_forces()
>>> forces.shape
(2, 3)
>>> print(f'Force: {forces[0, 2]:.3f} eV/Å')
Force: -0.639 eV/Å

.. image:: https://badge.fury.io/py/gpaw.svg
    :target: https://pypi.org/project/gpaw/

.. _Python: https://www.python.org
.. _ASE: https://wiki.fysik.dtu.dk/ase


.. _news:

News
====

* :ref:`GPAW version 24.6.0 <releasenotes>` released (Jun 7, 2024).

* `Psi-k highlight of the month
  <https://psi-k.net/download/highlights/Highlight_157.pdf>`__ (Apr 1, 2024)

* New publication:
  :doi:`GPAW: An open Python package for electronic structure calculations
  <10.1063/5.0182685>` (Mar 7, 2024)

* :ref:`GPAW version 24.1.0 <releasenotes>` released (Jan 4, 2024).

* .. warning::

    **IMPORTANT**: A bug was found in PW-mode `\Gamma`-point only calculations.
    Please check :ref:`here <bug0>` if you have been afected by this.

* :ref:`GPAW version 23.9.1 <releasenotes>` released (Sep 15, 2023).

* :ref:`GPAW version 23.9.0 <releasenotes>` released (Sep 13, 2023).

* Monthly *response code* two-day sprints will start on the last Monday
  of the month and continue the next day (Aug 28. 2023).

* Monthly *general maintenance* one-day sprints will start on the Tuesday
  in the week after the monthly response sprints
  (this will typically be the first Tuesday of the month, but it
  can also be the second Tuesday) (Aug 28. 2023).

* :ref:`GPAW version 23.6.1 <releasenotes>` released (Jul 5, 2023).

* :ref:`GPAW version 23.6.0 <releasenotes>` released (Jun 9, 2023).

* :ref:`GPAW version 22.8.0 <releasenotes>` released (Aug 18, 2022).

* :ref:`GPAW version 22.1.0 <releasenotes>` released (Jan 12, 2022).

* :ref:`GPAW version 21.6.0 <releasenotes>` released (Jun 24, 2021).

* Slides from the "GPAW 2021 Users and developers meeting" are
  now available `here
  <https://www.cecam.org/workshop-details/1039#document_tab>`__
  (Jun 2, 2021).

* Upcoming workshop:  The
  `GPAW 2021 Users and developers meeting
  <https://www.cecam.org/workshop-details/1039>`__
  will be held online on June 1--4, 2021.
  See also announcement on `Psi-k
  <https://psi-k.net/events/gpaw-2021-users-and-developers-meeting-june-1-4/>`__
  (Mar 1, 2021).

* :ref:`GPAW version 21.1.0 <releasenotes>` released (Jan 18, 2021).

* :ref:`GPAW version 20.10.0 <releasenotes>` released (Oct 19, 2020).

* :ref:`GPAW version 20.1.0 <releasenotes>` released (Jan 30, 2020).

* :ref:`GPAW version 19.8.1 <releasenotes>` released (Aug 8, 2019).

* :ref:`GPAW version 19.8.0 <releasenotes>` released (Aug 1, 2019).

* :ref:`GPAW version 1.5.2 <releasenotes>` released (May 8, 2019).

* :ref:`GPAW version 1.5.1 <releasenotes>` released (Jan 23, 2019).

* :ref:`GPAW version 1.5.0 <releasenotes>` released (Jan 11, 2019).

* :ref:`GPAW version 1.4.0 <releasenotes>` released (May 29, 2018).

* :ref:`GPAW version 1.3.0 <releasenotes>` released (Oct 2, 2017).

* Supported by NOMAD_ (Mar 1, 2017)

  .. image:: static/NOMAD_Logo_supported_by.png
     :width: 100 px
     :target: NOMAD_

* Code-sprints moved to first Tuesday of every month (Feb 17, 2017)

* :ref:`GPAW version 1.2 <releasenotes>` released (Feb 7, 2017)

* It has been decided to have monthly GPAW/ASE code-sprints at DTU in Lyngby.
  The sprints will be the first Wednesday of every month starting December 7,
  2016 (Nov 11, 2016)

* Slides from the talks at :ref:`workshop16` are now available (Sep 5, 2016)

* :ref:`GPAW version 1.1 <releasenotes>` released (Jun 22, 2016)

* :ref:`GPAW version 1.0 <releasenotes>` released (Mar 18, 2016)

* Web-page now use the `Read the Docs Sphinx Theme
  <https://github.com/snide/sphinx_rtd_theme>`_ (Mar 18, 2016)

* :ref:`GPAW version 0.11 <releasenotes>` released (Jul 22, 2015)

* :ref:`GPAW version 0.10 <releasenotes>` released (Apr 8, 2014)

* GPAW is part of the `PRACE Unified European Application Benchmark Suite`_
  (Oct 17, 2013)

* May 21-23, 2013: :ref:`GPAW workshop <workshop>` at the Technical
  University of Denmark (Feb 8, 2013)

* Prof. Häkkinen has received `18 million CPU hour grant`_ for GPAW based
  research project (Nov 20, 2012)

* A new :ref:`setups` bundle released (Oct 26, 2012)

* :ref:`GPAW version 0.9 <releasenotes>` released (March 7, 2012)

* :ref:`GPAW version 0.8 <releasenotes>` released (May 25, 2011)

* GPAW is part of benchmark suite for `CSC's supercomputer procurement`_
  (Apr 19, 2011)

* New features: Calculation of the linear :ref:`dielectric response
  <df_theory>` of an extended system (RPA and ALDA kernels) and
  calculation of :ref:`rpa` (Mar 18, 2011)

* Massively parallel GPAW calculations presented at `PyCon 2011`_.
  See William Scullin's talk here: `Python for High Performance
  Computing`_ (Mar 12, 2011)

* :ref:`GPAW version 0.7.2 <releasenotes>` released (Aug 13, 2010)

* :ref:`GPAW version 0.7 <releasenotes>` released (Apr 23, 2010)

* GPAW is `\Psi_k` `scientific highlight of the month`_ (Apr 3, 2010)

* A third GPAW code sprint was successfully hosted at CAMD (Oct 20, 2009)

* :ref:`GPAW version 0.6 <releasenotes>` released (Oct 9, 2009)

* `QuantumWise <https://www.quantumwise.com>`_ adds GPAW-support to
  `Virtual NanoLab`_ (Sep 8, 2009)

* Join the new IRC channel ``#gpaw`` on FreeNode (Jul 15, 2009)

* :ref:`GPAW version 0.5 <releasenotes>` released (Apr 1, 2009)

* A new :ref:`setups` bundle released (Mar 27, 2009)

* A second GPAW code sprint was successfully hosted at CAMD (Mar 20, 2009)

* :ref:`GPAW version 0.4 <releasenotes>` released (Nov 13, 2008)

* The :ref:`tutorialsexercises` are finally ready for use in the `CAMd summer
  school 2008`_ (Aug 15, 2008)

* This site is now powered by Sphinx_ (Jul 31, 2008)

* GPAW is now based on numpy_ instead of of Numeric (Jan 22, 2008)

* :ref:`GPAW version 0.3 <releasenotes>` released (Dec 19, 2007)

* CSC_ is organizing a `GPAW course`_: "Electronic structure
  calculations with GPAW" (Dec 11, 2007)

* The `code sprint 2007`_ was successfully finished (Nov 16, 2007)

* The source code is now in the hands of SVN and Trac (Oct 22, 2007)

* A GPAW Sprint will be held on November 16 in Lyngby (Oct 18, 2007)

* Work on atomic basis-sets begun (Sep 25, 2007)

.. _numpy: https://numpy.scipy.org/
.. _CSC: https://www.csc.fi
.. _GPAW course: https://www.csc.fi/english/csc/courses/archive/gpaw-2008-01
.. _Sphinx: https://www.sphinx-doc.org
.. _CAMd summer school 2008: https://www.camd.dtu.dk/English/Events/CAMD_Summer_School_2008/Programme.aspx
.. _code sprint 2007: https://www.dtu.dk/Nyheder/Nyt_fra_Institutterne.aspx?guid={38B92D63-FB09-4DFA-A074-504146A2D678}
.. _Virtual NanoLab: https://www.quantumwise.com/products/12-products/28-atk-se-200906#GPAW
.. _scientific highlight of the month: https://www.psi-k.org/newsletters/News_98/Highlight_98.pdf
.. _pycon 2011: https://us.pycon.org/2011/schedule/presentations/226/
.. _Python for High Performance Computing: https://pycon.blip.tv/file/4881240/
.. _CSC's supercomputer procurement: https://www.csc.fi/english/pages/hpc2011
.. _18 million CPU hour grant: https://www.prace-ri.eu/PRACE-5thRegular-Call
.. _PRACE Unified European Application Benchmark Suite: https://www.prace-ri.eu/ueabs
.. _NOMAD: https://repository.nomad-coe.eu/


.. toctree::

   algorithms
   install
   documentation/documentation
   tutorialsexercises/tutorialsexercises
   setups/setups
   releasenotes
   contact
   faq
   devel/devel
   summerschools/summerschools
   workshops/workshops
   bugs
