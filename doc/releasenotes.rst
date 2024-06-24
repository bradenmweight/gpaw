.. _releasenotes:

=============
Release notes
=============


Git master branch
=================

:git:`master <>`.

* ...


Version 24.6.0
==============

May 31, 2024: :git:`24.6.0 <../24.6.0>`

* Minimum version requirements: Python 3.8, ASE 3.23.0.

* Added a constraint on our Numpy dependency: ``numpy<2``.  GPAW almost
  works with numpy-2, but not quite.

* Fixed a bug in the electron-phonon module which constructed the PAW
  corrections to the Hamiltonian with incorrect off-diagonal elements. Our
  testing has not revealed any critical changes in results; however, if you've
  been using the electron-phonon module, we still recommend that you verify
  your results.

* New 14 electron Cr PAW potential added to our :ref:`setup releases`.
  For high accuracy, it is recommented over the old 6-electron version
  (which is still the default).  You can use it by
  specifying ``setups={'Cr': '14'}`` (see also :ref:`manual_setups`).
  It has been generated with the following command::

    $ gpaw dataset Cr -sw -r2.0 -P3s,4s,3p,4p,3d,d,F -fPBE -t 14 -b

  There is also an LDA version of the potential.

* A new ``adjust_cell()`` function is added to replace the
  ``Cluster.minimal_box()``, and the Cluster object will be
  discontinued. The new  ``adjust_cell()`` function expands the
  ``minimal_box()``
  function but considers the ``pbc`` for periodic calculations and
  only adjusts the cell in the non-periodic directions.


.. _bug0:

Version 24.1.0
==============

Jan 4, 2024: :git:`24.1.0 <../24.1.0>`

.. warning::

   PW-mode `\Gamma`-point calculations could sometimes find a fake
   eigenstate with eigenvalue equal to exactly 0 eV.  With a plane-wave
   expansion of the real-valued wave functions:

   .. math::

      \sum_\mathbf{G} c_\mathbf{G} e^{i\mathbf{G}\cdot\mathbf{r}}

   we must have `c_\mathbf{G}=c_\mathbf{-G}^*`, so the `c_\mathbf{G}`
   coefficient should not have an imaginary
   part for `\mathbf{G}=(0,0,0)`.  This was violated when
   the initial guess of the wave functions came from random numbers.

   For systems with vacuum in the cell, the 0 eV state would be unoccupied,
   but for fully periodic systems, the 0 eV state could be occupied.
   So, if you have done

   * PW-mode calculations
   * fully periodic systems
   * only `\Gamma`-point sampling
   * ``nbands`` set to a number large enough to trigger random
     wave functions

   then please check your results to see if you are affected by this bug.

   Fixed in :mr:`2114`.

* Minimum version requirements: Python 3.8, ASE 3.23.0b1, NumPy 1.17.0,
  SciPy 1.6.0

* Functionality has been added to calculate various local properties of the
  magnetic sites of a crystal, see :ref:`sites`.

* Python 3.7 no longer supported.

* Calculations of ground and excited states in FD and PW modes can now be
  done using direct orbital optimization (see :ref:`directmin`). Use
  ``eigensolver='etdm-fdpw'`` (or import ``gpaw.directmin.etdm_fdpw.FDPWETDM``
  ). The LCAO implementation of direct optimization available also in previous
  versions can be used by specifying ``eigensolver='etdm-lcao'`` (or by
  importing ``gpaw.directmin.etdm_lcao.LCAOETDM``).

  Direct optimization can also be used to perform calculations with
  Perdew-Zunger self-interaction correction (PZ-SIC) for both FD, PW and LCAO
  modes. Use (FD and PW modes)::

      from gpaw import GPAW
      from gpaw.directmin.etdm_fdpw import FDPWETDM
      calc = GPAW(eigensolver=LCAOETDM(localizationtype='PM_PZ',
                                       functional={'name': 'PZ-SIC'}),
                  ...)

  or (LCAO mode)::

      from gpaw import GPAW
      from gpaw.directmin.etdm_lcao import LCAOETDM
      calc = GPAW(eigensolver=LCAOETDM(localizationtype='PM_PZ',
                                       functional={'name': 'PZ-SIC'}),
                  ...)

  For excited state calculations, use direct optimization together with
  :ref:`MOM <mom>` (available for all modes) or with
  :ref:`generalized mode following <do-gmf>` (available only for LCAO).

* A bug in spin polarized (ferromagnetic) GW+BSE calculations was fixed:
  :issue:`828`.

* A bug resulting in slight inaccuracies when calculating inner products
  between radial partial waves was fixed. This bug affected the calculation
  of spin magnetic moments inside PAW spheres and the Hubbard correction when
  it was applied to p-states. See :issue:`1068`.


Version 23.9.1
==============

Sep 15, 2023: :git:`23.9.1 <../23.9.1>`

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.17.0,
  SciPy 1.6.0

* Include new GPU ``.c``, ``.h`` and ``.cpp`` files in :git:`MANIFEST.in`
  (:issue:`975`).


Version 23.9.0
==============

Sep 13, 2023: :git:`23.9.0 <../23.9.0>`

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.17.0,
  SciPy 1.6.0

* In the future, it will become an error to not specify a
  :ref:`mode <manual_mode>` parameter for a DFT calculation.
  For now, users will get a warning when finite-difference mode is
  implicitly chosen.  Please change your scripts to avoid this error/warning.

* Removed the utility function: ``gpaw.utilities.ibz2bz.ibz2bz``.

* :class:`~gpaw.tddft.TDDFT` and :class:`~gpaw.lcaotddft.LCAOTDDFT` will
  now throw an error if the ground state contains point group symmetries

* We are now using Pytest-cache for our :func:`gpaw.test.conftest.gpw_files`
  fixture.

* New
  :meth:`~gpaw.new.ase_interface.ASECalculator.get_orbital_magnetic_moments`
  method: calculates the orbital magnetic moment vector for each atom.

* New experimental density mixer: ``MixerFull``.


Version 23.6.1
==============

Jul 5, 2023: :git:`23.6.1 <../23.6.1>`

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.17.0,
  SciPy 1.6.0

* Avoid deprecation warnings from Numpy-1.25 (:mr:`1770`, :mr:`1771`).


Version 23.6.0
==============

Jun 9, 2023: :git:`23.6.0 <../23.6.0>`

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.17.0,
  SciPy 1.6.0

* New :meth:`~gpaw.calculator.GPAW.new()` method for creating new ``GPAW``
  objects with tweaked input parameters.  This can often be used to replace
  the use of the ``set`` method which we are planning to remove.

* A bug was found (now fixed) in the :ref:`zfs` module.  Please redo
  calculations done with versions 22.1 and 22.8.

* A bug in the implementation of MGGA functionals was found: :issue:`674`.
  The kinetic-energy density was calculated from the irreducible part of
  the Brillouin zone, but it was not symmetrized as it should be.  This
  has now been fixed.

  .. warning::

     If you have done any MGGA calculations taking advantage of symmetries
     in order to reduce number of **k**-points then you should redo those
     calculations.  Sorry!

* Two other bugs affecting MGGA calculations were found (in the fix_density
  and diagonalize_full_hamiltonian routines), which are fixed by
  `!1417 <https://gitlab.com/gpaw/gpaw/-/merge_requests/1417>`_.

  .. warning::

     MGGA calculations using fix_density and/or diagonalize_full_hamiltonian
     should be rerun with these fixes.

* The stress tensor was implemented for MGGA functionals, and
  parallelization of MGGAs for large systems was improved.

* Local orbitals added in LCAO mode to construct effective
  tight-binding Hamiltonians: :ref:`los in lcao`, :ref:`los tutorial`.

* Missing factor of `2\pi` now included in RPA shift current:
  :func:`gpaw.nlopt.shift.get_shift`.

* Updated RPA-energy tutorial: :ref:`c2cu rpa`.

* New tutorial: :ref:`abinitiomd`.

* Added relative tolerance for force convergence. This is useful for geometry
  optimizations to adaptively converge forces. See :ref:`custom_convergence`.

* Experimental support for PW-mode calculations using a GPU: :ref:`gpu`.

* One can now specify the total energy convergence criterium in eV instead
  of eV / valence electron:
  ``convergence={'energy': Energy(tol=..., relative=False)}``.
  See the :class:`gpaw.convergence_criteria.Energy` class.

* The PW-mode now includes an ``interpolation`` flag.  See
  :class:`gpaw.wavefunctions.pw.PW` for details.

* The LCAO implementation of direct optimization for variational calculations
  of excited electronic states now includes
  :ref:`constrained optimization <mom>`
  useful for challenging charge transfer excited states.

* The :ref:`direct optimization generalized mode following method <do-gmf>`
  (DO-GMF) for variational calculations of excited electronic states was added.

* Updated electron-phonon coupling and Raman implementations and
  documentation. See :ref:`elphtheory` and :ref:`elphraman`.

  .. warning::

     Bugs in previous versions could have led to wrong relative intensities.
     Please regenerate the Raman tensor.


Version 22.8.0
==============

Aug 18, 2022: :git:`22.8.0 <../22.8.0>`

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.15.0,
  SciPy 1.2.0

* Updated :ref:`WSL installation instructions <wsl>`.

* New feature for the :ref:`gpaw symmetry <cli>` command:  Will show number of
  **k**-points in the IBZ.

* New :class:`~gpaw.convergence_criteria.MaxIter` convergence criterium:
  ``convergence={'maximum iterations': 200}``.  This will let a calculation
  converge after 200 steps unless it already converged before that.  This is
  useful for structure optimizations that start far from the minimum.

* New common interface to the implementation of both linear and nonlinear
  frequency grids in the response code, now passed as a single input to e.g.
  Chi0, DielectricFunction and G0W0. Explained in the :ref:`frequency grid`
  tutorial.

* :ref:`spinspiral calculations`.  See also
  :git:`~gpaw/test/spinspiral/test_h_chain.py`.

* :ref:`soc`.

* The GW0 feature has been removed.

* :ref:`LrTDDFT <lrtddft>` works now also with LCAO-mode wave functions.

* GLLBSC functional uses now automatically Fermi level as the reference
  energy (GLLBSCM behavior) when the system has no band gap.
  This resolves "GLLBSC error: HOMO is higher than LUMO" observed in some
  systems during SCF iterations.
  See `!854 <https://gitlab.com/gpaw/gpaw/-/merge_requests/854>`_ for details.

* Functionality to compute magnon dispersions for ferromagnets in the
  classical isotropic Heisenberg model has been added, see
  ``gpaw.response.heisenberg``

* A new module ``gpaw.response.mft``, see :ref:`mft`, has been added for the
  calculation of isotropic Heisenberg exchange parameters within a linear
  response formulation of the magnetic force theorem. The module depends on a
  novel ``SiteKernels`` interface, see ``gpaw.response.site_kernels``, to
  discretize the DFT description into magnetic sublattices.


Version 22.1.0
==============

Jan 12, 2022: :git:`22.1.0 <../22.1.0>`

.. important::

   This release contains some important bug-fixes:

   * Spin-polarized GW-calculations:  The bug was introduced in
     version 20.10.0 and also present in versions 21.1.0 and 21.6.0.

   * Bug in non self-consistent eigenvalues for hybrid functionals
     and spin-polarized systems.

   * Erroneous Hirshfeld-effective volumes for non-orthogonal cells.

   * Fix for latest numpy-1.22.0.

* Minimum version requirements: Python 3.7, ASE 3.22.1, NumPy 1.15.0,
  SciPy 1.2.0

* Python 3.7 or later is required now.

* One can now apply Hund's rule (``hund=True``) to systems containing
  more than one atom.  This is useful for finding ferro-magnetic states
  and often works better that using ``magmoms=[1, 1, ...]`` for the
  initial magnetic moments.

* :ref:`polarizability` tutorial.

* Variational calculations of molecules and periodic systems in LCAO mode can
  now be done using the :ref:`exponential transformation direct minimization
  (ETDM) <directmin>`::

      from gpaw import GPAW
      calc = GPAW(eigensolver='etdm',
                  occupations={'name': 'fixed-uniform'},
                  mixer={'backend': 'no-mixing'},
                  nbands='nao',
                  ...)

  The use of ETDM is particularly recommended in
  excited-state calculations using MOM (see :ref:`mom`).

* Constant magnetic field calculations can now be done:
  See :class:`gpaw.bfield.BField` and this example:
  :git:`gpaw/test/ext_potential/test_b_field.py`.

* :ref:`raman` calculations for extended systems using electron-phonon
  coupling are now implemented in the LCAO mode.

  * An example can be found under :ref:`elphraman`.

  * The electron-phonon code has been updated. It can now be avoided to load
    the whole supercell matrix into memory.

  * A routine to calculate dipole and nabla (momentum) matrix elements for
    LCAO wave functions has been added: :git:`gpaw/lcao/dipoletransition.py`

* You can now change all sorts of things about how the SCF cycle decides it
  is converged. You can specify new, non-default convergence keywords like
  ``work function`` or ``minimum iterations``, you can change how default
  convergence keywords behave (like changing how many past energies the
  ``energy`` criterion examines), and you can even write your own custom
  convergence criteria. See :ref:`custom_convergence`.

* The SCF output table has been simplified, and a letter "c" now appears
  next to converged items.

* Charged molecule calculations with PW-mode have been improved.  The
  Poisson equation is now solved in a way so that monopole interactions
  between cells correctly vanish.

* The hyperfine tensor CLI-tool no longer divides by total magnetic moment:
  :ref:`hyperfine`.

* The solvated jellium method (:class:`~gpaw.solvation.sjm.SJM`)---for
  constant-potential calculations in simulating
  electrochemical/electrified interfaces---has been thoroughly
  updated, and more thorough :ref:`documentation<sjm>` and
  :ref:`tutorials<solvated_jellium_method>` are now available. Al keywords
  now enter the :class:`~gpaw.solvation.sjm.SJM` calculator through the
  :literal:`sj` dictionary.

* Radiative emission (lifetimes, ...) are obtainable from
  real-time LCAO-TDDFT via the radiation-reaction potential.
  See the tutorial: :ref:`radiation_reaction_rttddft`.


Version 21.6.0
==============

Jun 24, 2021: :git:`21.6.0 <../21.6.0>`

* Corresponding ASE release: ASE-3.22.0.

* :ref:`resonant_raman_water` tutorial added.

* The :ref:`time-propagation TDDFT (fd-mode) <timepropagation>` calculator
  refactored and observer support generalized.

  * The dipole moment output and restart file parameters are
    deprecated; use the corresponding observers instead.
    See the updated :ref:`documentation <timepropagation>`.

  * The observers for :ref:`inducedfield` need now to be defined before
    the kick instead of after it.

  * Corresponding updates for :ref:`qsfdtd` and :ref:`hybridscheme`.

* It is now possible to calculate electronic circular dichroism spectra
  with real-time time-propagation TDDFT.
  See the tutorial: :ref:`circular_dichroism_rtddft`.

* The documentation and tutorial for :ref:`lrtddft2` updated.

* True occupation numbers are now printed in the text output for the
  Kohn–Sham states.  Previously, the printed occupation numbers were
  scaled by **k**-point weight.

* Calculations of excited states can now be performed with the :ref:`Maximum
  Overlap Method (MOM) <mom>`. Since calculations using MOM are variational,
  they provide atomic forces and can be used for excited-state geometry
  optimization and molecular dynamics.

* The Davidson eigensolver now uses ScaLAPACK for the
  `(2 N_{\text{bands}}) \times (2 N_{\text{bands}})` diagonalization step
  when ``parallel={'sl_auto':True}`` is used.

* Removed several old command-line options:
  ``--memory-estimate-depth``, ``--domain-decomposition``,
  ``--state-parallelization``, ``--augment-grids``,
  ``--buffer-size``, ``--profile``, ``--gpaw``, ``--benchmark-imports``.
  See :ref:`manual_parallel` and :ref:`profiling` for alternatives.
  Instead of ``--gpaw=df_dry_run=N``, use the ``--dry-run=N`` option
  (see :ref:`command line options`).

* Added documentation for :ref:`elph` and added support for
  spin-polarized systems.

* Implemented multiple orbital Hubbard U corrections (EX: for correction
  of both p and d orbitals on transition metals)

* There used to be two versions of the GPAW web-page which was quite
  confusing.  The https://wiki.fysik.dtu.dk/gpaw/dev/ web-page has now been
  dropped.  There is now only https://wiki.fysik.dtu.dk/gpaw/ and it documents
  the use of the in development version of GPAW.

* ``gpaw sbatch`` will now detect an active virtual environment (venv)
  and activate it in the job script.


Version 21.1.0
===============

Jan 18, 2021: :git:`21.1.0 <../21.1.0>`

* Corresponding ASE release: ASE-3.21.0.

* We now use GPAW's own (faster) implementation for LDA, PBE, revPBE, RPBE
  and PW91.  For most calculation the speedup is unimportant, but for our
  test-suites it gives a nice boost.  There can be small meV changes compared
  to the LibXC implementation.  If you want to use LibXC then use::

      from gpaw.xc.gga import GGA
      from gpaw.xc.libxc import LibXC
      calc = GPAW(xc=GGA(LibXC('PBE')), ...)

* New :ref:`zfs` module.

* New :ref:`scissors operator`.

* Nonlinear optical responses can now be calculated in the independent
  particle approximations. See the :ref:`nlo_tutorial` tutorial for how
  to use it to compute the second-harmonic generation and shift current
  spectra.

* New method for interpolating pseudo density to fine grids:
  :meth:`gpaw.utilities.ps2ae.PS2AE.get_pseudo_density`
  (useful for Bader analysis and other things).

* Now with contribution from "frozen" core: :ref:`hyperfine`.

* Change in parameters of :ref:`linear response TDDFT <lrtddft>`

* Improved relaxation in the excited states in parallel,
  see  :ref:`linear response TDDFT <lrtddft>`

* We now have a code coverage report updated every night.

* Plane-wave mode implementation of hybrid functionals can now be selected
  via a *dict*: ``xc={'name': ..., 'backend': 'pw'}``, where then name must be
  one of EXX, PBE0, HSE03, HSE06 or B3LYP.  The EXX fraction and damping
  parameter can also be given in the dict.


Version 20.10.0
===============

Oct 19, 2020: :git:`20.10.0 <../20.10.0>`

* Corresponding ASE release: ASE-3.20.1.

* New :func:`gpaw.spinorbit.soc_eigenstates` function.  Handles
  parallelization and uses symmetry.  Angles are given in degrees
  (was radians before).

* The ``gpaw.spinorbit.get_anisotropy()`` method has been removed.  Use the
  :func:`~gpaw.spinorbit.soc_eigenstates` function combined with the
  :meth:`~gpaw.spinorbit.BZWaveFunctions.calculate_band_energy` method.
  See this tutorial: :ref:`magnetic anisotropy`.

* Improvements on GLLBSC and other GLLB-type exchange-correlation potentials:

  * `Fix for periodic metallic systems
    <https://gitlab.com/gpaw/gpaw/-/merge_requests/651>`_

  * `General fixes and improvements
    <https://gitlab.com/gpaw/gpaw/-/merge_requests/700>`_.
    Syntax for the discontinuity and band gap calculations has also been
    updated. See :ref:`the updated tutorial <band_gap>` for a detailed
    description of these calculations.

* Forces are now available for hybrid functionals in
  plane-wave mode.

* New functions for non self-consistent hybrid calculations:
  :func:`gpaw.hybrids.energy.non_self_consistent_energy` and
  :func:`gpaw.hybrids.eigenvalues.non_self_consistent_eigenvalues`.

* Python 3.6 or later is required now.

* Updates in :ref:`LCAOTDDFT <lcaotddft>` module:

  * User-defined time-dependent potentials and general kicks supported.

  * New observers for analysis.

  * Syntax updates for Kohn--Sham decomposition,
    see :ref:`examples <ksdecomposition>`.

  * Code improvements.

* New :meth:`~gpaw.calculator.GPAW.get_atomic_electrostatic_potentials`
  method.  Useful for aligning eigenvalues from different calculations.
  See :ref:`this example <potential>`.

* We are using pytest_ for testing.  Read about special GPAW-fixtures here:
  :ref:`testing`.

* We are now using MyPy_ for static analysis of the source code.

* Parallelization over spin is no longer possible.  This simplifies
  the code for handling non-collinear spins and spin-orbit coupling.

* Code for calculating occupation numbers has been refactored.  New functions:
  :func:`~gpaw.occupations.fermi_dirac`,
  :func:`~gpaw.occupations.marzari_vanderbilt` and
  :func:`~gpaw.occupations.methfessel_paxton`.  Deprecated:
  :func:`~gpaw.occupations.occupation_numbers`.  See :ref:`smearing`
  and :ref:`manual_occ` for details.

* Calculations with fixed occupation numbers are now done with
  ``occupations={'name': 'fixed', 'numbers': ...}``.

* The ``fixdensity`` keyword has been deprecated.

* New :meth:`gpaw.calculator.GPAW.fixed_density` method added to replace use
  of the deprecated ``fixdensity`` keyword.

* New configuration option (``nolibxc = True``) for compiling GPAW
  without LibXC.  This is mostly for debugging.  Only functionals available
  are: LDA, PBE, revPBE, RPBE and PW91.

* Tetrahedron method for Brillouin-zone integrations (**experimental**).
  Use ``occupations={'name': 'tetrahedron-method'}`` or
  ``occupations={'name': 'improved-tetrahedron-method'}``.
  See :doi:`Blöchl et. al <10.1103/PhysRevB.49.16223>`
  and :ref:`smearing` for details.

* New :func:`gpaw.mpi.broadcast_array` function for broadcasting
  an ``np.ndarray`` across several MPI-communicators.  New
  :func:`gpaw.mpi.send` and :func:`gpaw.mpi.receive` functions for general
  Python objects.

* Atoms with fractional atomic numbers can now be handled.

* When creating a ``GPAW`` calculator object from a gpw-file, the ``txt``
  defaults to ``None``.  Use ``GPAW('abc.gpw', txt='-')`` to get the old
  behavior.

* :ref:`hyperfine`.

* New :mod:`gpaw.point_groups` module.  See this tutorial:
  :ref:`point groups`.

* Default mixer (see :ref:`densitymix`) for spin-polarized systems has been
  changed from ``MixerSum`` to ``MixerDif``.  Now, both the total density
  and the magnetization density are mixed compared to before where only
  the total density was mixed.  To get the
  old behavior, use ``mixer=MixerSum(beta=0.05, history=5, weight=50)``
  for periodic systems
  and ``mixer=MixerSum(beta=0.25, history=3, weight=1)`` for molecules.

* New :func:`~gpaw.utilities.dipole.dipole_matrix_elements` and
  :func:`~gpaw.utilities.dipole.dipole_matrix_elements_from_calc`
  functions.  Command-line interface::

      $ python3 -m gpaw.utilities.dipole <gpw-file>


.. _pytest: https://doc.pytest.org/en/latest/contents.html
.. _mypy: https://mypy.readthedocs.io/en/stable/


Version 20.1.0
==============

Jan 30, 2020: :git:`20.1.0 <../20.1.0>`

* Corresponding ASE release: ASE-3.19.0.

* Self-consistent calculations with hybrid functionals are now possible in
  plane-wave mode.  You have to parallelize over plane-waves and you must
  use the Davidson eigensolver with one iteration per SCF step::

      from gpaw import GPAW, PW, Davidson
      calc = GPAW(mode=PW(ecut=...),
                  xc='HSE06',
                  parallel={'band': 1, 'kpt': 1},
                  eigensolver=Davidson(niter=1),
                  ...)

* We are now using setuptools_ instead of :mod:`distutils`.
  This means that installation with pip works much better.

* No more ``gpaw-python``.
  By default, an MPI-enabled Python interpreter is not built
  (use ``parallel_python_interpreter=True`` if you want a ``gpaw-python``).
  The ``_gpaw.so`` C-extension file (usually only used for serial calculations)
  will now be compiled with ``mpicc`` and contain what is necessary for both
  serial and parallel calculations.  In order to run GPAW in parallel, you
  do one of these three::

      $ mpiexec -n 24 gpaw python script.py
      $ gpaw -P 24 python script.py
      $ mpiexec -n 24 python3 script.py

  The first two are the recommended ones:  The *gpaw* script will make sure
  that imports are done in an efficient way.

* Configuration/customization:
  The ``customize.py`` file in the root folder of the Git repository is no
  longer used.  Instead, the first of the following three files that exist
  will be used:

  1) the file that ``$GPAW_CONFIG`` points at
  2) ``<git-root>/siteconfig.py``
  3) ``~/.gpaw/siteconfig.py``

  This will be used to configure things
  (BLAS, FFTW, ScaLAPACK, libxc, libvdwxc, ...).  If no configuration file
  is found then you get ``libraries = ['xc', 'blas']``.

* A Lapack library is no longer needed for compiling GPAW.  We are using
  :mod:`scipy.linalg` from now on.

* Debug mode is now enabled with::

      $ python3 -d script.py

* Dry-run mode is now enabled with::

      $ gpaw python --dry-run=N script.py

* New convergence criterium.  Example: ``convergence={'bands': 'CBM+2.5'}``
  will converge bands up to conduction band minimum plus 2.5 eV.

* Point-group symmetries now also used for non-periodic systems.
  Use ``symmetry={'point_group': False}`` if you don't want that.

* :ref:`Marzari-Vanderbilt distribution function <manual_occ>` added.

* New configuration option: ``noblas = True``.  Useful for compiling GPAW
  without a BLAS library.  :mod:`scipy.linalg.blas` and :func:`numpy.dot`
  will be used instead.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/


Version 19.8.1
==============

Aug 8, 2019: :git:`19.8.1 <../19.8.1>`

.. warning:: Upgrading from version 1.5.2

    Some small changes in the code introduced between version 1.5.2 and
    19.8.1 (improved handling of splines) may give rise to small changes in
    the total energy calculated with version 19.8.1 compared
    to version 1.5.2.  The changes should be in the meV/atom range, but may
    add up to significant numbers if you are doing calculations for large
    systems with many atoms.

* Corresponding ASE release: ASE-3.18.0.

* *Important bug fixed*: reading of some old gpw-files did not work.


Version 19.8.0
==============

Aug 1, 2019: :git:`19.8.0 <../19.8.0>`

* Corresponding ASE release: ASE-3.18.0.

* The ``"You have a weird unit cell"`` and
  ``"Real space grid not compatible with symmetry operation"``
  errors are now gone.  GPAW now handles these cases by
  choosing the number of real-space grid-points in a more clever way.

* The angular part of the PAW correction to the ALDA kernel is now calculated
  analytically by expanding the correction in spherical harmonics.

* Berry phases can now be calculated.  See the :ref:`berry tutorial` tutorial
  for how to use it to calculate spontaneous polarization, Born effective
  charges and other physical properties.

* How to do :ref:`ehrenfest` has now been documented.

* Non self-consistent hybrid functional calculations can now be continued if
  they run out of time.

* When using a convergence criteria on the accuracy of the forces
  (see :ref:`manual_convergence`), the forces will only be calculated when the
  other convergence criteria (energy, eigenstates and density) are fulfilled.
  This can save a bit of time.

* Experimental support for JTH_ PAW-datasets.

* Fast C implementation of bond-length constraints and associated hidden
  constraints for water models. This allows efficient explicit solvent QMMM
  calculations for GPAW up to tens of thousands of solvent molecules with
  water models such as SPC, TIPnP etc.  See :git:`gpaw/utilities/watermodel.py`
  and :git:`gpaw/test/test_rattle.py` for examples.

* New "metallic boundary conditions" have been added to the for PoissonSolver.
  This enables simulating charged 2D systems without counter charges.
  See: :git:`gpaw/test/poisson/test_metallic_poisson.py`

* Removed unnecessary application of H-operator in Davidson algorithm making
  it a bit faster.

.. _JTH: https://www.abinit.org/psp-tables


Version 1.5.2
=============

May 8, 2019: :git:`1.5.2 <../1.5.2>`

* Corresponding ASE release: ASE-3.17.0.

* **Important bugfix release**:

  There was a bug which was triggered when combining
  ScaLAPACK, LCAO and k-points in GPAW 1.5.0/1.5.1 from January.  The
  projections were calculated incorrectly which affected the SCF
  loop.

  If you use ScaLAPACK+LCAO+kpoints and see the line "Atomic Correction:
  distributed and sparse using scipy" in the output, then please rerun
  after updating.


Version 1.5.1
=============

Jan 23, 2019: :git:`1.5.1 <../1.5.1>`

* Corresponding ASE release: ASE-3.17.0.

* Small bug fixes related to latest versions of Python, Numpy and Libxc.


Version 1.5.0
=============

Jan 11, 2019: :git:`1.5.0 <../1.5.0>`

* Corresponding ASE release: ASE-3.17.0.

* Last release to support Python 2.7.

* The default finite-difference stencils used for gradients in GGA and MGGA
  calculations have been changed.

  * The range of the stencil has been increased
    from 1 to 2 thereby decreasing the error from `O(h^2)` to `O(h^4)`
    (where `h` is the grid spacing).  Use ``xc={'name': 'PBE', 'stencil': 1}``
    to get the old, less accurate, stencil.

  * The stencils are now symmetric also for non-orthorhombic
    unit cells.  Before, the stencils would only have weight on the
    neighboring grid-points in the 6 directions along the lattice vectors.
    Now, grid-points along all nearest neighbor directions can have a weight
    in the  stencils.  This allows for creating stencils that have all the
    crystal symmetries.

* PW-mode calculations can now be parallelized over plane-wave coefficients.

* The PW-mode code is now much faster.  The "hot spots" have been moved
  from Python to C-code.

* Wavefunctions are now updated when the atomic positions change by
  default, improving the initial wavefunctions across geometry steps.
  Corresponds to ``GPAW(experimental={'reuse_wfs_method': 'paw'})``.
  To get the old behavior, set the option to ``'keep'`` instead.
  The option is disabled for TDDFT/Ehrenfest.

* Add interface to ELPA eigensolver for LCAO mode.
  Using ELPA is strongly recommended for large calculations.
  Use::

      GPAW(mode='lcao',
           basis='dzp',
           parallel={'sl_auto': True, 'use_elpa': True})

  See also documentation on the :ref:`parallel keyword <manual_parallel>`.

* Default eigensolver is now ``Davidson(niter=2)``.

* Default number of bands is now `1.2 \times N_{\text{occ}} + 4`, where
  `N_{\text{occ}}` is the number of occupied bands.

* Solvated jellium method has been implemented, see
  :ref:`the documentation <solvated_jellium_method>`.

* Added FastPoissonSolver which is faster and works well for any cell.
  This replaces the old Poisson solver as default Poisson solver.

* :ref:`rsf` and improved virtual orbitals, the latter from Hartree-Fock
  theory.

* New Jupyter notebooks added for teaching DFT and many-body methods.  Topics
  cover: :ref:`catalysis`, :ref:`magnetism`, :ref:`machinelearning`,
  :ref:`photovoltaics`, :ref:`batteries` and :ref:`intro`.

* New experimental local **k**-point refinement feature:
  :git:`gpaw/test/test_kpt_refine.py`.

* A module and tutorial have been added for calculating electrostatic
  corrections to DFT total energies for charged systems involving localized
  defects: :ref:`defects`.

* Default for FFTW planning has been changed from ``ESTIMATE`` to ``MEASURE``.
  See :class:`gpaw.wavefunctions.pw.PW`.


Version 1.4.0
=============

May 29, 2018: :git:`1.4.0 <../1.4.0>`

* Corresponding ASE release: ASE-3.16.0.

* Improved parallelization of operations with localized functions in
  PW mode.  This solves the current size bottleneck in PW mode.

* Added QNA XC functional: :ref:`qna`.

* Major refactoring of the LCAOTDDFT code and added Kohn--Sham decomposition
  analysis within LCAOTDDFT, see :ref:`the documentation <lcaotddft>`.

* New ``experimental`` keyword, ``GPAW(experimental={...})`` to enable
  features that are still being tested.

* Experimental support for calculations with non-collinear spins
  (plane-wave mode only).
  Use ``GPAW(experimental={'magmoms': magmoms})``, where ``magmoms``
  is an array of magnetic moment vectors of shape ``(len(atoms), 3)``.

* Number of bands no longer needs to be divisible by band parallelization
  group size.  Number of bands will no longer be automatically adjusted
  to fit parallelization.

* Major code refactoring to facilitate work with parallel arrays.  See new
  module: ``gpaw.matrix``.

* Better reuse of wavefunctions when atoms are displaced.  This can
  improve performance of optimizations and dynamics in FD and PW mode.
  Use ``GPAW(experimental={'reuse_wfs_method': name})`` where name is
  ``'paw'`` or ``'lcao'``.  This will move the projections of the
  wavefunctions upon the PAW projectors or LCAO basis set along with
  the atoms.  The latter is best when used with ``dzp``.
  This feature has no effect for LCAO mode where the basis functions
  automatically follow the atoms.

* Broadcast imports (Python3 only): Master process broadcasts most module
  files at import time to reduce file system overhead in parallel
  calculations.

* Command-line arguments for BLACS/ScaLAPACK
  have been
  removed in favor of the :ref:`parallel keyword
  <manual_parallelization_types>`.  For example instead of running
  ``gpaw-python --sl_diagonalize=4,4,64``, set the parallelization
  within the script using
  ``GPAW(parallel={'sl_diagonalize': (4, 4, 64)})``.

* When run through the ordinary Python interpreter, GPAW will now only
  intercept and use command-line options of the form ``--gpaw
  key1=value1,key2=value2,...`` or ``--gpaw=key1=value1,key2=value2,...``.

* ``gpaw-python`` now takes :ref:`command line options` directly
  instead of stealing them from ``sys.argv``, passing the remaining
  ones to the script:
  Example: ``gpaw-python --gpaw=debug=True myscript.py myscript_arguments``.
  See also ``gpaw-python --help``.

* Two new parameters for specifying the Pulay stress. Directly like this::

      GPAW(mode=PW(ecut, pulay_stress=...), ...)

  or indirectly::

      GPAW(mode=PW(ecut, dedecut=...), ...)

  via the formula `\sigma_P=(2/3)E_{\text{cut}}dE/dE_{\text{cut}}/V`.  Use
  ``dedecut='estimate'`` to use an estimate from the kinetic energy of an
  isolated atom.

* New utility function: ``gpaw.utilities.ibz2bz.ibz2bz``.


Version 1.3.0
=============

October 2, 2017: :git:`1.3.0 <../1.3.0>`

* Corresponding ASE release: ASE-3.15.0.

* :ref:`command line options` ``--dry-run`` and ``--debug`` have been removed.
  Please use ``--gpaw dry-run=N`` and ``--gpaw debug=True`` instead
  (or ``--gpaw dry-run=N,debug=True`` for both).

* The :meth:`ase.Atoms.get_magnetic_moments` method will no longer be
  scaled to sum up to the total magnetic moment.  Instead, the magnetic
  moments integrated inside the atomic PAW spheres will be returned.

* New *sbatch* sub-command for GPAW's :ref:`cli`.

* Support added for ASE's new *band-structure* :ref:`ase:cli`::

  $ ase band-structure xxx.gpw -p GKLM

* Added :ref:`tetrahedron method <tetrahedron>` for calculation the density
  response function.

* Long-range cutoff for :mod:`~ase.calculators.qmmm` calculations can now be
  per molecule instead of only per point charge.

* Python 2.6 no longer supported.

* There is now a web-page documenting the use of the in development version
  of GPAW: https://wiki.fysik.dtu.dk/gpaw/dev/.

* :ref:`BSE <bse tutorial>` calculations for spin-polarized systems.

* Calculation of :ref:`magnetic anisotropy <magnetic anisotropy>`.

* Calculation of vectorial magnetic moments inside PAW spheres based on
  spin-orbit spinors.

* Added a simple :func:`gpaw.occupations.occupation_numbers` function for
  calculating occupation numbers, Fermi-level, magnetic moment, and entropy
  from eigenvalues and k-point weights.

* Deprecated calculator-keyword ``dtype``.  If you need to force the datatype
  of the wave functions to be complex, then use something like::

      calc = GPAW(mode=PW(ecut=500, force_complex_dtype=True))

* Norm-conserving potentials (HGH and SG15) now subtract the Hartree
  energies of the compensation charges.
  The total energy of an isolated pseudo-atom stripped of all valence electrons
  will now be zero.

* HGH and SG15 pseudopotentials are now Fourier-filtered at run-time
  as appropriate for the given grid spacing.  Using them now requires scipy.

* The ``gpaw dos`` sub-command of the :ref:`cli` can now show projected DOS.
  Also, one can now use linear tetrahedron interpolation for the calculation
  of the (P)DOS.

* The :class:`gpaw.utilities.ps2ae.PS2AE` tool can now also calculate the
  all-electron electrostatic potential.


Version 1.2.0
=============

Feb 7, 2017: :git:`1.2.0 <../1.2.0>`.

* Corresponding ASE release: ASE-3.13.0.

* New file-format for gpw-files.  Reading of old files should still work.
  Look inside the new files with::

      $ python3 -m ase.io.ulm abc.gpw

* Simple syntax for specifying BZ paths introduced:
  ``kpts={'path': 'GXK', 'npoints': 50}``.

* Calculations with ``fixdensity=True`` no longer update the Fermi level.

* The GPAW calculator object has a new
  :meth:`~ase.calculators.calculator.Calculator.band_structure`
  method that returns an :class:`ase.spectrum.band_structure.BandStructure`
  object.  This makes it easy to create band-structure plots as shown
  in section 9 of this awesome Psi-k *Scientfic Highlight Of The Month*:
  https://psi-k.net/download/highlights/Highlight_134.pdf.

* Dipole-layer corrections for slab calculations can now be done in PW-mode
  also.  See :ref:`dipole`.

* New :meth:`~gpaw.calculator.GPAW.get_electrostatic_potential` method.

* When setting the default PAW-datasets or basis-sets using a dict, we
  must now use ``'default'`` as the key instead of ``None``:

  >>> calc = GPAW(basis={'default': 'dzp', 'H': 'sz(dzp)'})

  and not:

  >>> calc = GPAW(basis={None: 'dzp', 'H': 'sz(dzp)'})

  (will still work, but you will get a warning).

* New feature added to the GW code to be used with 2D systems. This lowers
  the required k-point grid necessary for convergence. See this tutorial
  :ref:`gw-2D`.

* It is now possible to carry out GW calculations with eigenvalue self-
  consistency in G. NOTE: This feature was removed after version 22.1.0.

* XC objects can now be specified as dictionaries, allowing GGAs and MGGAs
  with custom stencils: ``GPAW(xc={'name': 'PBE', 'stencil': 2})``

* Support for spin-polarized vdW-DF functionals (svdW-DF) with libvdwxc.


Version 1.1.0
=============

June 22, 2016: :git:`1.1.0 <../1.1.0>`.

* Corresponding ASE release: ASE-3.11.0.

* There was a **BUG** in the recently added spin-orbit module.  Should now
  be fixed.

* The default Davidson eigensolver can now parallelize over bands.

* There is a new PAW-dataset file available:
  :ref:`gpaw-setup-0.9.20000.tar.gz <datasets>`.
  It's identical to the previous
  one except for one new data-file which is needed for doing vdW-DF
  calculations with Python 3.

* Jellium calculations can now be done in plane-wave mode and there is a new
  ``background_charge`` keyword (see the :ref:`Jellium tutorial <jellium>`).

* New band structure unfolding tool and :ref:`tutorial <unfolding tutorial>`.

* The :meth:`~gpaw.calculator.GPAW.get_pseudo_wave_function` method
  has a new keyword:  Use ``periodic=True`` to get the periodic part of the
  wave function.

* New tool for interpolating the pseudo wave functions to a fine real-space
  grids and for adding PAW-corrections in order to obtain all-electron wave
  functions.  See this tutorial: :ref:`ps2ae`.

* New and improved dataset pages (see :ref:`periodic table`).  Now shows
  convergence of absolute and relative energies with respect to plane-wave
  cut-off.

* :ref:`wannier90 interface`.

* Updated MacOSX installation guide for :ref:`homebrew` users.

* topological index


Version 1.0.0
=============

Mar 17, 2016: :git:`1.0.0 <../1.0.0>`.

* Corresponding ASE release: ASE-3.10.0.

* A **BUG** related to use of time-reversal symmetry was found in the
  `G_0W_0` code that was introduced in version 0.11.  This has been `fixed
  now`_ --- *please run your calculations again*.

* New :mod:`gpaw.external` module.

* The gradients of the cavity and the dielectric in the continuum
  solvent model are now calculated analytically for the case of the
  effective potential method. This improves the accuracy of the forces
  in solution compared to the gradient calculated by finite
  differences. The solvation energies are expected to change slightly
  within the accuracy of the model.

* New `f_{\text{xc}}` kernels for correlation energy calculations.  See this
  updated :ref:`tutorial <rapbe_tut>`.

* Correlation energies within the range-separated RPA.

* Experimental interface to the libvdwxc_ library
  for efficient van der Waals density functionals.

* It's now possible to use Davidson and CG eigensolvers for MGGA calculations.

* The functional name "M06L" is now deprecated.  Use "M06-L" from now on.


.. _fixed now: https://gitlab.com/gpaw/gpaw/commit/c72e02cd789
.. _libvdwxc: https://gitlab.com/libvdwxc/libvdwxc


Version 0.11.0
==============

July 22, 2015: :git:`0.11.0 <../0.11.0>`.

* Corresponding ASE release: ASE-3.9.1.

* When searching for basis sets, the setup name if any is now
  prepended automatically to the basis name.  Thus if
  :file:`setups='<setupname>'` and :file:`basis='<basisname>'`, GPAW
  will search for :file:`<symbol>.<setupname>.<basisname>.basis`.

* :ref:`Time-propagation TDDFT with LCAO <lcaotddft>`.

* Improved distribution and load balance when calculating atomic XC
  corrections, and in LCAO when calculating atomic corrections to the
  Hamiltonian and overlap.

* Norm-conserving :ref:`SG15 pseudopotentials <manual_setups>` and
  parser for several dialects of the UPF format.

* Non self-consistent spin-orbit coupling have been added. See :ref:`tutorial
  <spinorbit>` for examples of band structure calculations with spin-orbit
  coupling.

* Text output from ground-state calculations now list the symmetries found
  and the **k**-points used.  Eigenvalues and occupation numbers are now
  also printed for systems with **k**-points.

* :ref:`GW <gw exercise>`, :ref:`rpa`, and :ref:`response function
  calculation <df_tutorial>` has been rewritten to take advantage of
  symmetry and fast matrix-matrix multiplication (BLAS).

* New :ref:`symmetry <manual_symmetry>` keyword.  Replaces ``usesymm``.

* Use non-symmorphic symmetries: combining fractional translations with
  rotations, reflections and inversion.  Use
  ``symmetry={'symmorphic': False}`` to turn this feature on.

* New :ref:`forces <manual_convergence>` keyword in convergence.  Can
  be used to calculate forces to a given precision.

* Fixed bug in printing work functions for calculations with a
  dipole-correction `<https://listserv.fysik.dtu.dk/pipermail/
  gpaw-users/2015-February/003226.html>`_.

* A :ref:`continuum solvent model <continuum_solvent_model>` was added.

* A :ref:`orbital-free DFT <ofdft>` with PAW transformation is available.

* GPAW can now perform :ref:`electrodynamics` simulations using the
  quasistatic finite-difference time-domain (QSFDTD) method.

* BEEF-vdW, mBEEF and mBEEF-vdW functionals added.

* Support for Python 3.


Version 0.10.0
==============

Apr 8, 2014: :git:`0.10.0 <../0.10.0>`.

* Corresponding ASE release: ASE-3.8.1

* Default eigensolver is now the Davidson solver.

* Default density mixer parameters have been changed for calculations
  with periodic boundary conditions.  Parameters for that case:
  ``Mixer(0.05, 5, 50)`` (or ``MixerSum(0.05, 5, 50)`` for spin-paired
  calculations).  Old parameters: ``0.1, 3, 50``.

* Default is now ``occupations=FermiDirac(0.1)`` if a
  calculation is periodic in at least one direction,
  and ``FermiDirac(0.0)`` otherwise (before it was 0.1 eV for anything
  with **k**-points, and 0 otherwise).

* Calculations with a plane-wave basis set are now officially supported.

* :ref:`One-shot GW calculations <gw_theory>` with full frequency
  integration or plasmon-pole approximation.

* Beyond RPA-correlation: `using renormalized LDA and PBE
  <https://trac.fysik.dtu.dk/projects/gpaw/browser/branches/sprint2013/doc/tutorials/fxc_correlation>`_.

* :ref:`bse theory`.

* Improved RMM-DIIS eigensolver.

* Support for new libxc 2.0.1.  libxc must now be built separately from GPAW.

* MGGA calculations can be done in plane-wave mode.

* Calculation of the stress tensor has been implemented for plane-wave
  based calculation (except MGGA).

* MGGA: number of neighbor grid points to use for FD stencil for
  wave function gradient changed from 1 to 3.

* New setups: Y, Sb, Xe, Hf, Re, Hg, Tl, Rn

* Non self-consistent calculations with screened hybrid functionals
  (HSE03 and HSE06) can be done in plane-wave mode.

* Modified setups:

  .. note::

     Most of the new semi-core setups currently require
     :ref:`eigensolver <manual_eigensolver>` ``dav``, ``cg``
     eigensolvers or ``rmm-diis`` eigensolver with a couple of iterations.

  - improved egg-box: N, O, K, S, Ca, Sc, Zn, Sr, Zr, Cd, In, Sn, Pb, Bi

  - semi-core states included: Na, Mg, V, Mn, Ni,
    Nb, Mo, Ru (seems to solve the Ru problem :git:`gpaw/test/big/Ru001/`),
    Rh, Pd, Ag, Ta, W, Os, Ir, Pt

  - semi-core states removed: Te

  - elements removed: La (energetics was wrong: errors ~1eV per unit cell
    for PBE formation energy of La2O3 wrt. PBE benchmark results)

  .. note::

     For some of the setups one has now a choice of different
     number of valence electrons, e.g.::

       setups={'Ag': '11'}

     See :ref:`manual_setups` and list the contents of
     :envvar:`GPAW_SETUP_PATH` for available setups.

* new ``dzp`` basis set generated for all the new setups, see
  https://trac.fysik.dtu.dk/projects/gpaw/ticket/241


Version 0.9.0
=============

Mar 7, 2012: :git:`0.9.0 <../0.9.0>`.

* Corresponding ASE release: ase-3.6

* Convergence criteria for eigenstates changed: The missing volume per
  grid-point factor is now included and the units are now eV**2. The
  new default value is 4.0e-8 eV**2 which is equivalent to the old
  default for a grid spacing of 0.2 Å.

* GPAW should now work also with NumPy 1.6.

* Much improved :ref:`cli` now based on the `new tool`_ in ASE.


.. _new tool: https://wiki.fysik.dtu.dk/ase/ase/cmdline.html


Version 0.8.0
=============

May 25, 2011: :git:`0.8.0 <../0.8.0>`.

* Corresponding ASE release: ase-3.5.1
* Energy convergence criterion changed from 1 meV/atom to 0.5
  meV/electron.  This was changed in order to allow having no atoms like
  for jellium calculations.
* Linear :ref:`dielectric response <df_theory>` of an extended system
  (RPA and ALDA kernels) can now be calculated.
* :ref:`rpa`.
* Non self-consistent calculations with k-points for hybrid functionals.
* Methfessel-Paxton distribution added.
* Text output now shows the distance between planes of grid-points as
  this is what will be close to the grid-spacing parameter *h* also for
  non-orthorhombic cells.
* Exchange-correlation code restructured.  Naming convention for
  explicitly specifying libxc functionals has changed: :ref:`manual_xc`.
* New PAW setups for Rb, Ti, Ba, La, Sr, K, Sc, Ca, Zr and Cs.


Version 0.7.2
=============

Aug 13, 2010: :git:`0.7.2 <../0.7.2>`.

* Corresponding ASE release: ase-3.4.1
* For version 0.7, the default Poisson solver was changed to
  ``PoissonSolver(nn=3)``.  Now, also the Poisson solver's default
  value for ``nn`` has been changed from ``'M'`` to ``3``.


Version 0.7
===========

Apr 23, 2010: :git:`0.7 <../0.7>`.

* Corresponding ASE release: ase-3.4.0
* Better and much more efficient handling of non-orthorhombic unit
  cells.  It may actually work now!
* Much better use of ScaLAPACK and BLACS.  All large matrices can now
  be distributed.
* New test coverage pages for all files.
* New default value for Poisson solver stencil: ``PoissonSolver(nn=3)``.
* Much improved MPI module (:ref:`communicators`).
* Self-consistent Meta GGA.
* New :ref:`PAW setup tar-file <setups>` now contains revPBE setups and
  also dzp basis functions.
* New ``$HOME/.gpaw/rc.py`` configuration file.
* License is now GPLv3+.
* New HDF IO-format.
* :ref:`Advanced GPAW Test System <big-test>` Introduced.


Version 0.6
===========

Oct 9, 2009: :git:`0.6 <../0.6>`.

* Corresponding ASE release: ase-3.2.0
* Much improved default parameters.
* Using higher order finite-difference stencil for kinetic energy.
* Many many other improvements like: better parallelization, fewer bugs and
  smaller memory footprint.


Version 0.5
===========

Apr 1, 2009: :git:`0.5 <../0.5>`.

* Corresponding ASE release: ase-3.1.0
* `new setups added Bi, Br, I, In, Os, Sc, Te; changed Rb setup <https://trac.fysik.dtu.dk/projects/gpaw/changeset/3612>`_.
* `memory estimate feature is back <https://trac.fysik.dtu.dk/projects/gpaw/changeset/3575>`_


Version 0.4
===========

Nov 13, 2008: :git:`0.4 <../0.4>`.

* Corresponding ASE release: ase-3.0.0
* Now using ASE-3 and numpy.
* TPSS non self-consistent implementation.
* LCAO mode.
* vdW-functional now coded in C.
* Added atomic orbital basis generation scripts.
* Added an Overlap object, and moved ``apply_overlap`` and
  ``apply_hamiltonian`` from ``Kpoint`` to Overlap and Hamiltonian classes.

* Wannier code much improved.
* Experimental LDA+U code added.
* Now using libxc.
* Many more setups.
* Delta SCF calculations.

* Using localized functions will now no longer use MPI group
  communicators and blocking calls to MPI_Reduce and MPI_Bcast.
  Instead non-blocking sends/receives/waits are used.  This will
  reduce synchronization time for large parallel calculations.
* More work on LB94.
* Using LCAO code for initial guess for grid calculations.
* TDDFT.
* Moved documentation to Sphinx.
* Improved metric for Pulay mixing.
* Porting and optimization for BlueGene/P.
* Experimental Hartwigsen-Goedecker-Hutter pseudopotentials added.
* Transport calculations with LCAO.


Version 0.3
===========

Dec 19, 2007: :git:`0.3 <../0.3>`.
