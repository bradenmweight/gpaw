.. _faq:

==========================
Frequently Asked Questions
==========================

.. contents::


General
=======

.. _citation:

Citation: How should I cite GPAW?
---------------------------------

If you find GPAW useful in your research please cite this GPAW review:

   | Jens Jørgen Mortensen, Ask Hjorth Larsen, Mikael Kuisma et al.
   | :doi:`GPAW: An open Python package for electronic structure calculations featured <10.1063/5.0182685>`
   | J. Chem. Phys. 160, 092503 (2024)

together with the ASE review (see :ref:`ase:cite`).

You are welcome to cite also the original GPAW reference and an earlier GPAW review:

   | J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen
   | :doi:`Real-space grid implementation of the projector augmented wave method <10.1103/PhysRevB.71.035109>`
   | Phys. Rev. B **71**, 035109 (2005)

   | J. Enkovaara, C. Rostgaard, J. J. Mortensen et al.
   | :doi:`Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method <10.1088/0953-8984/22/25/253202>`
   | J. Phys.: Condens. Matter **22**, 253202 (2010)

Please also cite those of the following that are relevant to you work:

* `Libxc <https://libxc.gitlab.io>`_ for XC-functionals other
  than LDA, PBE, revPBE, RPBE and PW91:

    *S. Lehtola, C. Steigemann, M. J. T. Oliveira and M. A. L. Marques.*,
    :doi:`Recent developments in LIBXC — a comprehensive library of functionals
    for density functional theory <10.1016/j.softx.2017.11.002>`,
    SoftwareX **7**, 1 (2018)

* :ref:`timepropagation` or :ref:`lrtddft`:

    *M. Walter, H. Häkkinen, L. Lehtovaara, M. Puska, J. Enkovaara,
    C. Rostgaard and J. J. Mortensen*,
    :doi:`Time-dependent density-functional theory in the
    projector augmented-wave method <10.1063/1.2943138>`,
    J. Chem. Phys. **128**, 244101 (2008)

* :ref:`Localized basis set calculations <lcao>` (LCAO):

    *A. H. Larsen, M. Vanin, J. J. Mortensen, K. S. Thygesen, and
    K. W. Jacobsen*,
    :doi:`Localized atomic basis set in the projector augmented wave
    method <10.1103/PhysRevB.80.195112>`,
    Phys. Rev. B **80**, 195112 (2009)

* :ref:`Linear dielectric response of an extended systems <df_tutorial>`:

    *J. Yan, J. J. Mortensen, K. W. Jacobsen, and K. S. Thygesen*,
    :doi:`Linear density response function in the
    projector augmented wave method: Applications to solids, surfaces,
    and interfaces <10.1103/PhysRevB.83.245122>`,
    Phys. Rev. B **83**, 245122 (2011)

* :ref:`Quasi-particle spectrum in the GW approximation <gw tutorial>`:

    *F. Hüser, T. Olsen, and K. S. Thygesen*,
    :doi:`Quasiparticle GW calculations for solids, molecules,
    and two-dimensional materials <10.1103/PhysRevB.87.235132>`,
    Phys. Rev. B **87**, 235132 (2013)

* :ref:`continuum_solvent_model`:

    *A. Held and M. Walter*,
    :doi:`Simplified continuum solvent model with a smooth cavity based
    on volumetric data <10.1063/1.4900838>`,
    J. Chem. Phys. **141**, 174108 (2014)

* :ref:`lcaotddft`:

    *M. Kuisma, A. Sakko, T. P. Rossi, A. H. Larsen, J. Enkovaara,
    L. Lehtovaara, and T. T. Rantala*,
    :doi:`Localized surface plasmon resonance in silver nanoparticles:
    Atomistic first-principles time-dependent density functional theory
    calculations <10.1103/PhysRevB.91.115431>`,
    Phys. Rev. B **91**, 115431 (2015)

* :ref:`ksdecomposition` and :ref:`lcaotddft`:

    *T. P. Rossi, M. Kuisma, M. J. Puska, R. M. Nieminen, and P. Erhart*,
    :doi:`Kohn--Sham Decomposition in Real-Time Time-Dependent
    Density-Functional Theory: An Efficient Tool for Analyzing
    Plasmonic Excitations <10.1021/acs.jctc.7b00589>`,
    J. Chem. Theory Comput. **13**, 4779 (2017)

* :ref:`sjm`:

    *G. Kastlunger, P. Lindgren, A.A. Peterson*,
    :doi:`Controlled-potential simulation of elementary electrochemical
    reactions: proton discharge on metal surfaces <10.1021/acs.jpcc.8b02465>`,
    J. Phys. Chem. C **122**, 12771 (2018)


Citations of the GPAW method papers
-----------------------------------

.. image:: documentation/citations.png
   :width: 750

(updated on 18 Mar 2021)

The total number of citations above is the number of publications
citing at least one of the other papers, not the sum of all citation
counts.

BibTex (:git:`doc/GPAW.bib`):

.. literalinclude:: GPAW.bib
   :language: bibtex


How do you pronounce GPAW?
--------------------------

In English: "geepaw" with a long "a".

In Danish: Først bogstavet "g", derefter "pav": "g-pav".

In Finnish: supisuomalaisittain "kee-pav".

In Polish: "gyeh" jak `"Gie"rek <https://en.wikipedia.org/wiki/Edward_Gierek>`_, "pav" jak `paw <https://pl.wikipedia.org/wiki/Paw_indyjski>`_: "gyeh-pav".


Compiling the C-code
====================

For architecture dependent settings see the
:ref:`platforms and architectures` page.

Compilation of the C part failed::

 [~]$ python2.4 setup.py build_ext
 building '_gpaw' extension
 pgcc -fno-strict-aliasing -DNDEBUG -O2 -g -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -m64 -D_GNU_SOURCE -fPIC -fPIC -I/usr/include/python2.4 -c c/localized_functions.c -o build/temp.linux-x86_64-2.4/c/localized_functions.o -Wall -std=c99
 pgcc-Warning-Unknown switch: -fno-strict-aliasing
 PGC-S-0040-Illegal use of symbol, _Complex (/usr/include/bits/cmathcalls.h: 54)

You are probably using another compiler, than was used for
compiling python. Undefine the environment variables CC, CFLAGS and
LDFLAGS with::

 # sh/bash users:
 unset CC; unset CFLAGS; unset LDFLAGS
 # csh/tcsh users:
 unsetenv CC; unsetenv CFLAGS; unsetenv LDFLAGS

and try again.


Calculation does not converge
=============================

Consult the :ref:`convergence` page.


Poisson solver did not converge!
================================

If you are doing a spin-polarized calculation for an isolated molecule,
then you should set the Fermi temperature to a low value.

You can also try to set the number of grid points to be divisible by 8.
Consult the :ref:`poisson_performance` page.
