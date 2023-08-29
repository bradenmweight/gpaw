import os
from contextlib import contextmanager
from pathlib import Path
import functools

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.lattice.hexagonal import Graphene
from ase.io import read
from gpaw import GPAW, PW, Davidson, FermiDirac, setup_paths
from gpaw.cli.info import info
from gpaw.mpi import broadcast, world
from gpaw.utilities import devnull
from ase.lattice.compounds import L1_2
from gpaw import Mixer
from gpaw.new.ase_interface import GPAW as GPAWNew


@contextmanager
def execute_in_tmp_path(request, tmp_path_factory):
    if world.rank == 0:
        # Obtain basename as
        # * request.function.__name__  for function fixture
        # * request.module.__name__    for module fixture
        basename = getattr(request, request.scope).__name__
        path = tmp_path_factory.mktemp(basename)
    else:
        path = None
    path = broadcast(path)
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield path
    finally:
        os.chdir(cwd)


@pytest.fixture(scope='function')
def in_tmp_dir(request, tmp_path_factory):
    """Run test function in a temporary directory."""
    with execute_in_tmp_path(request, tmp_path_factory) as path:
        yield path


@pytest.fixture(scope='module')
def module_tmp_path(request, tmp_path_factory):
    """Run test module in a temporary directory."""
    with execute_in_tmp_path(request, tmp_path_factory) as path:
        yield path


@pytest.fixture
def add_cwd_to_setup_paths():
    """Temporarily add current working directory to setup_paths."""
    try:
        setup_paths[:0] = ['.']
        yield
    finally:
        del setup_paths[:1]


response_band_cutoff = dict(
)


def with_band_cutoff(*, gpw, band_cutoff):
    # Store the band cutoffs in a dictionary to aid response tests
    response_band_cutoff[gpw] = band_cutoff

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, band_cutoff=band_cutoff, **kwargs)
        return wrapper

    return decorator


@pytest.fixture(scope='session')
def gpw_files(request):
    """Reuse gpw-files.

    Returns a dict mapping names to paths to gpw-files.
    The files are written to the pytest cache and can be cleared using
    pytest --cache-clear.

    Example::

        def test_something(gpw_files):
            calc = GPAW(gpw_files['h2_lcao'])
            ...

    Possible systems are:

    * Bulk BCC-Li with 3x3x3 k-points: ``bcc_li_pw``, ``bcc_li_fd``,
      ``bcc_li_lcao``.

    * O2 molecule: ``o2_pw``.

    * H2 molecule: ``h2_pw``, ``h2_fd``, ``h2_lcao``.

    * H2 molecule (not centered): ``h2_pw_0``.

    * Spin-polarized H atom: ``h_pw``.

    * Polyethylene chain.  One unit, 3 k-points, no symmetry:
      ``c2h4_pw_nosym``.  Three units: ``c6h12_pw``.

    * Bulk TiO2 with 4x4x4 k-points: ``ti2o4_pw`` and ``ti2o4_pw_nosym``.

    * Bulk BN (zinkblende) with 2x2x2 k-points and 9 converged bands:
      ``bn_pw``.

    * h-BN layer with 3x3x1 (gamma center) k-points and 26 converged bands:
      ``hbn_pw``.

    * Graphene with 6x6x1 k-points: ``graphene_pw``

    * MoS2 with 6x6x1 k-points: ``mos2_pw`` and ``mos2_pw_nosym``

    * NiCl2 with 6x6x1 k-points: ``nicl2_pw``

    * V2Br4 (AFM monolayer), LDA, 4x2x1 k-points, 28(+1) converged bands:
      ``v2br4_pw`` and ``v2br4_pw_nosym``

    * Bulk Si, LDA, 2x2x2 k-points (gamma centered): ``si_pw``

    * Bulk Si, LDA, 4x4x4 k-points, 8(+1) converged bands: ``fancy_si_pw``
      and ``fancy_si_pw_nosym``

    * Bulk Fe, LDA, 4x4x4 k-points, 9(+1) converged bands: ``fe_pw``
      and ``fe_pw_nosym``

    * Bulk C, LDA, 2x2x2 k-points (gamma centered), ``c_pw``

    * Bulk Co (HCP), 4x4x4 k-points, 12(+1) converged bands: ``co_pw``
      and ``co_pw_nosym``

    * Bulk SrVO3 (SC), 3x3x3 k-points, 20(+1) converged bands: ``srvo3_pw``
      and ``srvo3_pw_nosym``

    * Bulk Al, LDA, 4x4x4 k-points, 10(+1) converged bands: ``al_pw``
      and ``al_pw_nosym``

    * Bulk Al, LDA, 4x4x4 k-points, 4 converged bands: ``bse_al``

    * Bulk Ag, LDA, 2x2x2 k-points, 6 converged bands,
      2eV U on d-band: ``ag_pw``

    * Bulk GaAs, LDA, 4x4x4 k-points, 8(+1) bands converged: ``gaas_pw``
      and ``gaas_pw_nosym``

    * Distorted bulk Fe, revTPSS: ``fe_pw_distorted``

    * Distorted bulk Si, TPSS: ``si_pw_distorted``

    Files always include wave functions.
    """
    cache = request.config.cache
    gpaw_cachedir = cache.mkdir('gpaw_test_gpwfiles')
    return GPWFiles(gpaw_cachedir)


class Locked(FileExistsError):
    pass


@contextmanager
def temporary_lock(path):
    fd = None
    try:
        with path.open('x') as fd:
            yield
    except FileExistsError:
        raise Locked()
    finally:
        if fd is not None:
            path.unlink()


@contextmanager
def world_temporary_lock(path):
    if world.rank == 0:
        try:
            with temporary_lock(path):
                world.sum_scalar(1)
                yield
        except Locked:
            world.sum_scalar(0)
            raise
    else:
        status = world.sum_scalar(0)
        if status:
            yield
        else:
            raise Locked


_all_gpw_methodnames = set()


def gpwfile(meth):
    """Decorator to identify the methods that produce gpw files."""
    _all_gpw_methodnames.add(meth.__name__)
    return meth


class GPWFiles:
    """Create gpw-files."""
    def __init__(self, path: Path):
        self.path = path

        self.gpw_files = {}
        for file in path.glob('*.gpw'):
            self.gpw_files[file.name[:-4]] = file

    def __getitem__(self, name: str) -> Path:
        if name in self.gpw_files:
            return self.gpw_files[name]

        gpwpath = self.path / (name + '.gpw')

        lockfile = self.path / f'{name}.lock'

        for _attempt in range(60):  # ~60s timeout
            files_exist = 0
            if world.rank == 0:
                files_exist = int(gpwpath.exists())
            files_exist = world.sum_scalar(files_exist)

            if files_exist:
                self.gpw_files[name] = gpwpath
                return self.gpw_files[name]

            try:
                with world_temporary_lock(lockfile):
                    calc = getattr(self, name)()
                    work_path = gpwpath.with_suffix('.tmp')
                    calc.write(work_path, mode='all')

                    # By now files should exist *and* be fully written, by us.
                    # Rename them to the final intended paths:
                    if world.rank == 0:
                        work_path.rename(gpwpath)

            except Locked:
                import time
                time.sleep(1)

        raise RuntimeError(f'GPW fixture generation takes too long: {name}.  '
                           'Consider using pytest --cache-clear if there are '
                           'stale lockfiles, else write faster tests.')

    @gpwfile
    def bcc_li_pw(self):
        return self.bcc_li({'name': 'pw', 'ecut': 200})

    @gpwfile
    def bcc_li_fd(self):
        return self.bcc_li({'name': 'fd'})

    @gpwfile
    def bcc_li_lcao(self):
        return self.bcc_li({'name': 'lcao'})

    def bcc_li(self, mode):
        li = bulk('Li', 'bcc', 3.49)
        li.calc = GPAW(mode=mode,
                       kpts=(3, 3, 3),
                       txt=self.path / f'bcc_li_{mode["name"]}.txt')
        li.get_potential_energy()
        return li.calc

    @gpwfile
    def fcc_Ni_col(self):
        return self.fcc_Ni('col')

    @gpwfile
    def fcc_Ni_ncol(self):
        return self.fcc_Ni('ncol')

    @gpwfile
    def fcc_Ni_ncolsoc(self):
        return self.fcc_Ni('ncolsoc')

    def fcc_Ni(self, calc_type):
        Ni = bulk('Ni', 'fcc', 3.48)
        Ni.center()

        mm = 0.5
        easy_axis = 1 / np.sqrt(3) * np.ones(3)
        Ni.set_initial_magnetic_moments([mm])

        symmetry = {'point_group': True, 'time_reversal': True} if \
            calc_type == 'col' else 'off'
        magmoms = None if calc_type == 'col' else [mm * easy_axis]
        soc = True if calc_type == 'ncolsoc' else False

        Ni.calc = GPAWNew(mode={'name': 'pw', 'ecut': 400}, xc='LDA',
                          kpts={'size': (4, 4, 4), 'gamma': True},
                          parallel={'domain': 1, 'band': 1},
                          symmetry=symmetry,
                          occupations={'name': 'fermi-dirac', 'width': 0.05},
                          magmoms=magmoms, soc=soc,
                          txt=self.path / f'fcc_Ni_{calc_type}.txt')
        Ni.get_potential_energy()
        return Ni.calc

    @gpwfile
    def h2_pw(self):
        return self.h2({'name': 'pw', 'ecut': 200})

    @gpwfile
    def h2_fd(self):
        return self.h2({'name': 'fd'})

    @gpwfile
    def h2_lcao(self):
        return self.h2({'name': 'lcao'})

    def h2(self, mode):
        h2 = Atoms('H2', positions=[[0, 0, 0], [0.74, 0, 0]])
        h2.center(vacuum=2.5)
        h2.calc = GPAW(mode=mode,
                       txt=self.path / f'h2_{mode["name"]}.txt')
        h2.get_potential_energy()
        return h2.calc

    @gpwfile
    def h2_pw_0(self):
        h2 = Atoms('H2',
                   positions=[[-0.37, 0, 0], [0.37, 0, 0]],
                   cell=[5.74, 5, 5],
                   pbc=True)
        h2.calc = GPAW(mode={'name': 'pw', 'ecut': 200},
                       txt=self.path / 'h2_pw_0.txt')
        h2.get_potential_energy()
        return h2.calc

    @gpwfile
    def h2_bcc_afm(self):
        a = 2.75
        atoms = bulk(name='H', crystalstructure='bcc', a=a, cubic=True)
        atoms.set_initial_magnetic_moments([1., -1.])

        atoms.calc = GPAW(xc='LDA',
                          txt=self.path / 'h2_bcc_afm.txt',
                          mode=PW(250),
                          nbands=4,
                          convergence={'bands': 4},
                          kpts={'density': 2.0, 'gamma': True})
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def h_pw(self):
        h = Atoms('H', magmoms=[1])
        h.center(vacuum=4.0)
        h.calc = GPAW(mode={'name': 'pw', 'ecut': 500},
                      txt=self.path / 'h_pw.txt')
        h.get_potential_energy()
        return h.calc

    @gpwfile
    def o2_pw(self):
        d = 1.1
        a = Atoms('O2', positions=[[0, 0, 0], [d, 0, 0]], magmoms=[1, 1])
        a.center(vacuum=4.0)
        a.calc = GPAW(mode={'name': 'pw', 'ecut': 800},
                      txt=self.path / 'o2_pw.txt')
        a.get_potential_energy()
        return a.calc

    @gpwfile
    def Cu3Au_qna(self):
        ecut = 300
        kpts = (1, 1, 1)

        QNA = {'alpha': 2.0,
               'name': 'QNA',
               'stencil': 1,
               'orbital_dependent': False,
               'parameters': {'Au': (0.125, 0.1), 'Cu': (0.0795, 0.005)},
               'setup_name': 'PBE',
               'type': 'qna-gga'}

        atoms = L1_2(['Au', 'Cu'], latticeconstant=3.7)
        atoms[0].position[0] += 0.01  # Break symmetry already here
        calc = GPAW(mode=PW(ecut),
                    eigensolver=Davidson(2),
                    nbands='120%',
                    mixer=Mixer(0.4, 7, 50.0),
                    parallel=dict(domain=1),
                    convergence={'density': 1e-4},
                    xc=QNA,
                    kpts=kpts,
                    txt=self.path / 'Cu3Au.txt')
        atoms.calc = calc
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def co_lcao(self):
        d = 1.1
        co = Atoms('CO', positions=[[0, 0, 0], [d, 0, 0]])
        co.center(vacuum=4.0)
        co.calc = GPAW(mode='lcao',
                       txt=self.path / 'co_lcao.txt')
        co.get_potential_energy()
        return co.calc

    @gpwfile
    def c2h4_pw_nosym(self):
        d = 1.54
        h = 1.1
        x = d * (2 / 3)**0.5
        z = d / 3**0.5
        pe = Atoms('C2H4',
                   positions=[[0, 0, 0],
                              [x, 0, z],
                              [0, -h * (2 / 3)**0.5, -h / 3**0.5],
                              [0, h * (2 / 3)**0.5, -h / 3**0.5],
                              [x, -h * (2 / 3)**0.5, z + h / 3**0.5],
                              [x, h * (2 / 3)**0.5, z + h / 3**0.5]],
                   cell=[2 * x, 0, 0],
                   pbc=(1, 0, 0))
        pe.center(vacuum=2.0, axis=(1, 2))
        pe.calc = GPAW(mode='pw',
                       kpts=(3, 1, 1),
                       symmetry='off',
                       txt=self.path / 'c2h4_pw_nosym.txt')
        pe.get_potential_energy()
        return pe.calc

    @gpwfile
    def c6h12_pw(self):
        pe = read(self['c2h4_pw_nosym'])
        pe = pe.repeat((3, 1, 1))
        pe.calc = GPAW(mode='pw', txt=self.path / 'c6h12_pw.txt')
        pe.get_potential_energy()
        return pe.calc

    @gpwfile
    def h2o_lcao(self):
        from ase.build import molecule
        atoms = molecule('H2O', cell=[8, 8, 8], pbc=1)
        atoms.center()
        atoms.calc = GPAW(mode='lcao', txt=self.path / 'h2o.txt')
        atoms.get_potential_energy()
        return atoms.calc

    def ti2o4(self, symmetry):
        pwcutoff = 400.0
        k = 4
        a = 4.59
        c = 2.96
        u = 0.305

        rutile_cell = [[a, 0, 0],
                       [0, a, 0],
                       [0, 0, c]]

        TiO2_basis = np.array([[0.0, 0.0, 0.0],
                               [0.5, 0.5, 0.5],
                               [u, u, 0.0],
                               [-u, -u, 0.0],
                               [0.5 + u, 0.5 - u, 0.5],
                               [0.5 - u, 0.5 + u, 0.5]])

        bulk_crystal = Atoms(symbols='Ti2O4',
                             scaled_positions=TiO2_basis,
                             cell=rutile_cell,
                             pbc=(1, 1, 1))

        tag = '_nosym' if symmetry == 'off' else ''
        bulk_calc = GPAW(mode=PW(pwcutoff),
                         nbands=42,
                         eigensolver=Davidson(1),
                         kpts={'size': (k, k, k), 'gamma': True},
                         xc='PBE',
                         occupations=FermiDirac(0.00001),
                         parallel={'band': 1},
                         symmetry=symmetry,
                         txt=self.path / f'ti2o4_pw{tag}.txt')

        bulk_crystal.calc = bulk_calc
        bulk_crystal.get_potential_energy()
        return bulk_calc

    @gpwfile
    def ti2o4_pw(self):
        return self.ti2o4({})

    @gpwfile
    def ti2o4_pw_nosym(self):
        return self.ti2o4('off')

    @gpwfile
    def si_pw(self):
        si = bulk('Si')
        calc = GPAW(mode='pw',
                    xc='LDA',
                    occupations=FermiDirac(width=0.001),
                    kpts={'size': (2, 2, 2), 'gamma': True},
                    txt=self.path / 'si_pw.txt')
        si.calc = calc
        si.get_potential_energy()
        return si.calc

    @gpwfile
    @with_band_cutoff(gpw='fancy_si_pw',
                      band_cutoff=8)  # 2 * (3s, 3p)
    def _fancy_si(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}
        xc = 'LDA'
        kpts = 4
        pw = 300
        occw = 0.01
        conv = {'bands': band_cutoff + 1,
                'density': 1.e-8}
        atoms = bulk('Si')
        atoms.center()

        tag = '_nosym' if symmetry == 'off' else ''
        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw),
            kpts={'size': (kpts, kpts, kpts), 'gamma': True},
            nbands=band_cutoff + 12,  # + 2 * (4s, 3d),
            occupations=FermiDirac(occw),
            convergence=conv,
            txt=self.path / f'fancy_si_pw{tag}.txt',
            symmetry=symmetry)

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def fancy_si_pw(self):
        return self._fancy_si()

    @gpwfile
    def fancy_si_pw_nosym(self):
        return self._fancy_si(symmetry='off')

    @gpwfile
    def bn_pw(self):
        atoms = bulk('BN', 'zincblende', a=3.615)
        atoms.calc = GPAW(mode=PW(400),
                          kpts={'size': (2, 2, 2), 'gamma': True},
                          nbands=12,
                          convergence={'bands': 9},
                          occupations=FermiDirac(0.001),
                          txt=self.path / 'bn_pw.txt')
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def hbn_pw(self):
        atoms = Graphene(symbol='B',
                         latticeconstant={'a': 2.5, 'c': 1.0},
                         size=(1, 1, 1))
        atoms[0].symbol = 'N'
        atoms.pbc = (1, 1, 0)
        atoms.center(axis=2, vacuum=3.0)
        atoms.calc = GPAW(mode=PW(400),
                          xc='LDA',
                          nbands=50,
                          occupations=FermiDirac(0.001),
                          parallel={'domain': 1},
                          convergence={'bands': 26},
                          kpts={'size': (3, 3, 1), 'gamma': True})
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def graphene_pw(self):
        from ase.lattice.hexagonal import Graphene
        atoms = Graphene(symbol='C',
                         latticeconstant={'a': 2.45, 'c': 1.0},
                         size=(1, 1, 1))
        atoms.pbc = (1, 1, 0)
        atoms.center(axis=2, vacuum=4.0)
        ecut = 250
        nkpts = 6
        atoms.calc = GPAW(mode=PW(ecut),
                          kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                          nbands=len(atoms) * 6,
                          txt=self.path / 'graphene_pw.txt')
        atoms.get_potential_energy()
        return atoms.calc

    def _mos2(self, symmetry=None):
        if symmetry is None:
            symmetry = {}
        from ase.build import mx2
        atoms = mx2(formula='MoS2', kind='2H', a=3.184, thickness=3.127,
                    size=(1, 1, 1), vacuum=5)
        atoms.pbc = (1, 1, 0)
        ecut = 250
        nkpts = 6
        tag = '_nosym' if symmetry == 'off' else ''
        atoms.calc = GPAW(mode=PW(ecut),
                          xc='LDA',
                          kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                          occupations=FermiDirac(0.01),
                          txt=self.path / f'mos2_pw{tag}.txt',
                          symmetry=symmetry)

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def mos2_pw(self):
        return self._mos2()

    @gpwfile
    def mos2_pw_nosym(self):
        return self._mos2(symmetry='off')

    @gpwfile
    def ni_pw_kpts333(self):
        from ase.dft.kpoints import monkhorst_pack
        # from gpaw.mpi import serial_comm
        Ni = bulk('Ni', 'fcc')
        Ni.set_initial_magnetic_moments([0.7])

        kpts = monkhorst_pack((3, 3, 3))

        calc = GPAW(mode='pw',
                    kpts=kpts,
                    occupations=FermiDirac(0.001),
                    setups={'Ni': '10'},
                    parallel=dict(domain=1),  # >1 fails on 8 cores
                    # communicator=serial_comm
                    )

        Ni.calc = calc
        Ni.get_potential_energy()
        calc.diagonalize_full_hamiltonian()
        # calc.write('Ni.gpw', mode='all')
        return calc

    @gpwfile
    def c_pw(self):
        atoms = bulk('C')
        atoms.center()
        calc = GPAW(mode=PW(150),
                    convergence={'bands': 6},
                    nbands=12,
                    kpts={'gamma': True, 'size': (2, 2, 2)},
                    xc='LDA')

        atoms.calc = calc
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def nicl2_pw(self):
        from ase.build import mx2

        # Define input parameters
        xc = 'LDA'
        kpts = 6
        pw = 300
        occw = 0.01
        conv = {'density': 1.e-8,
                'forces': 1.e-8}

        a = 3.502
        thickness = 2.617
        vacuum = 3.0
        mm = 2.0

        # Set up atoms
        atoms = mx2(formula='NiCl2', kind='1T', a=a,
                    thickness=thickness, vacuum=vacuum)
        atoms.set_initial_magnetic_moments([mm, 0.0, 0.0])
        # Use pbc to allow for real-space density interpolation
        atoms.pbc = True

        # Set up calculator
        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw,
                    # Interpolate the density in real-space
                    interpolation=3),
            kpts={'size': (kpts, kpts, 1), 'gamma': True},
            occupations=FermiDirac(occw),
            convergence=conv,
            txt=self.path / 'nicl2_pw.txt')

        atoms.get_potential_energy()

        return atoms.calc

    @with_band_cutoff(gpw='v2br4_pw',
                      band_cutoff=28)  # V(4s,3d) = 6, Br(4s,4p) = 4
    def _v2br4(self, *, band_cutoff, symmetry=None):
        from ase.build import mx2

        if symmetry is None:
            symmetry = {}

        # Define input parameters
        xc = 'LDA'
        kpts = 4
        pw = 200
        occw = 0.01
        conv = {'density': 1.e-4,
                'bands': band_cutoff + 1}

        a = 3.840
        thickness = 2.897
        vacuum = 3.0
        mm = 3.0

        # Set up atoms
        atoms = mx2(formula='VBr2', kind='1T', a=a,
                    thickness=thickness, vacuum=vacuum)
        atoms = atoms.repeat((1, 2, 1))
        atoms.set_initial_magnetic_moments([mm, 0.0, 0.0, -mm, 0.0, 0.0])
        # Use pbc to allow for real-space density interpolation
        atoms.pbc = True

        # Set up calculator
        tag = '_nosym' if symmetry == 'off' else ''
        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw,
                    # Interpolate the density in real-space
                    interpolation=3),
            kpts={'size': (kpts, kpts // 2, 1), 'gamma': True},
            setups={'V': '5'},
            nbands=band_cutoff + 12,
            occupations=FermiDirac(occw),
            convergence=conv,
            symmetry=symmetry,
            txt=self.path / f'v2br4_pw{tag}.txt')

        atoms.get_potential_energy()

        return atoms.calc

    @gpwfile
    def v2br4_pw(self):
        return self._v2br4()

    @gpwfile
    def v2br4_pw_nosym(self):
        return self._v2br4(symmetry='off')

    @with_band_cutoff(gpw='fe_pw',
                      band_cutoff=9)  # 4s, 4p, 3d = 9
    def _fe(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}
        """See also the fe_fixture_test.py test."""
        xc = 'LDA'
        kpts = 4
        pw = 300
        occw = 0.01
        conv = {'bands': band_cutoff + 1,
                'density': 1.e-8}
        a = 2.867
        mm = 2.21
        atoms = bulk('Fe', 'bcc', a=a)
        atoms.set_initial_magnetic_moments([mm])
        atoms.center()
        tag = '_nosym' if symmetry == 'off' else ''

        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw),
            kpts={'size': (kpts, kpts, kpts)},
            nbands=band_cutoff + 9,
            occupations=FermiDirac(occw),
            convergence=conv,
            txt=self.path / f'fe_pw{tag}.txt',
            symmetry=symmetry)

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def fe_pw(self):
        return self._fe()

    @gpwfile
    def fe_pw_nosym(self):
        return self._fe(symmetry='off')

    @with_band_cutoff(gpw='co_pw',
                      band_cutoff=14)  # 2 * (4s + 3d)
    def _co(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}
        # ---------- Inputs ---------- #

        # Atomic configuration
        a = 2.5071
        c = 4.0695
        mm = 1.6
        atoms = bulk('Co', 'hcp', a=a, c=c)
        atoms.set_initial_magnetic_moments([mm, mm])
        atoms.center()

        # Ground state parameters
        xc = 'LDA'
        occw = 0.01
        ebands = 2 * 2  # extra bands for ground state calculation
        pw = 200
        conv = {'density': 1e-8,
                'forces': 1e-8,
                'bands': band_cutoff + 1}

        # ---------- Calculation ---------- #

        tag = '_nosym' if symmetry == 'off' else ''
        atoms.calc = GPAW(xc=xc,
                          mode=PW(pw),
                          kpts={'size': (4, 4, 4), 'gamma': True},
                          occupations=FermiDirac(occw),
                          convergence=conv,
                          nbands=band_cutoff + ebands,
                          symmetry=symmetry,
                          txt=self.path / f'co_pw{tag}.txt')

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def co_pw(self):
        return self._co()

    @gpwfile
    def co_pw_nosym(self):
        return self._co(symmetry='off')

    @with_band_cutoff(gpw='srvo3_pw',
                      band_cutoff=20)
    def _srvo3(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}

        nk = 3
        cell = bulk('V', 'sc', a=3.901).cell
        atoms = Atoms('SrVO3', cell=cell, pbc=True,
                      scaled_positions=((0.5, 0.5, 0.5),
                                        (0, 0, 0),
                                        (0, 0.5, 0),
                                        (0, 0, 0.5),
                                        (0.5, 0, 0)))
        # Ground state parameters
        xc = 'LDA'
        occw = 0.01
        ebands = 10  # extra bands for ground state calculation
        pw = 200
        conv = {'density': 1e-8,
                'bands': band_cutoff + 1}

        # ---------- Calculation ---------- #

        tag = '_nosym' if symmetry == 'off' else ''
        atoms.calc = GPAW(xc=xc,
                          mode=PW(pw),
                          kpts={'size': (nk, nk, nk), 'gamma': True},
                          occupations=FermiDirac(occw),
                          convergence=conv,
                          nbands=band_cutoff + ebands,
                          symmetry=symmetry,
                          txt=self.path / f'srvo3_pw{tag}.txt')

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def srvo3_pw(self):
        return self._srvo3()

    @gpwfile
    def srvo3_pw_nosym(self):
        return self._srvo3(symmetry='off')

    @with_band_cutoff(gpw='al_pw',
                      band_cutoff=10)  # 3s, 3p, 4s, 3d
    def _al(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}
        xc = 'LDA'
        kpts = 4
        pw = 300
        occw = 0.01
        conv = {'bands': band_cutoff + 1,
                'density': 1.e-8}
        a = 4.043
        atoms = bulk('Al', 'fcc', a=a)
        atoms.center()
        tag = '_nosym' if symmetry == 'off' else ''

        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw),
            kpts={'size': (kpts, kpts, kpts), 'gamma': True},
            nbands=band_cutoff + 4,  # + 4p, 5s
            occupations=FermiDirac(occw),
            convergence=conv,
            txt=self.path / f'al_pw{tag}.txt',
            symmetry=symmetry)

        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def al_pw(self):
        return self._al()

    @gpwfile
    def al_pw_nosym(self):
        return self._al(symmetry='off')

    @gpwfile
    def bse_al(self):
        a = 4.043
        atoms = bulk('Al', 'fcc', a=a)
        calc = GPAW(mode='pw',
                    kpts={'size': (4, 4, 4), 'gamma': True},
                    xc='LDA',
                    nbands=4,
                    convergence={'bands': 'all'})

        atoms.calc = calc
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def ag_plusU_pw(self):
        xc = 'LDA'
        kpts = 2
        nbands = 6
        pw = 300
        occw = 0.01
        conv = {'bands': nbands,
                'density': 1e-12}
        a = 4.07
        atoms = bulk('Ag', 'fcc', a=a)
        atoms.center()

        atoms.calc = GPAW(
            xc=xc,
            mode=PW(pw),
            kpts={'size': (kpts, kpts, kpts), 'gamma': True},
            setups={'Ag': '11:d,2.0,0'},
            nbands=nbands,
            occupations=FermiDirac(occw),
            convergence=conv,
            parallel={'domain': 1},
            txt=self.path / 'ag_pw.txt')

        atoms.get_potential_energy()

        atoms.calc.diagonalize_full_hamiltonian()

        return atoms.calc

    @gpwfile
    def gaas_pw_nosym(self):
        return self._gaas(symmetry='off')

    @gpwfile
    def gaas_pw(self):
        return self._gaas()

    @with_band_cutoff(gpw='gaas_pw',
                      band_cutoff=8)
    def _gaas(self, *, band_cutoff, symmetry=None):
        if symmetry is None:
            symmetry = {}
        nk = 4
        cell = bulk('Ga', 'fcc', a=5.68).cell
        atoms = Atoms('GaAs', cell=cell, pbc=True,
                      scaled_positions=((0, 0, 0), (0.25, 0.25, 0.25)))
        tag = '_nosym' if symmetry == 'off' else ''
        conv = {'bands': band_cutoff + 1,
                'density': 1.e-8}

        calc = GPAW(mode=PW(400),
                    xc='LDA',
                    occupations=FermiDirac(width=0.01),
                    convergence=conv,
                    nbands=band_cutoff + 1,
                    kpts={'size': (nk, nk, nk), 'gamma': True},
                    txt=self.path / f'gs_GaAs{tag}.txt',
                    symmetry=symmetry)

        atoms.calc = calc
        atoms.get_potential_energy()
        return atoms.calc

    @gpwfile
    def h_pw280_fulldiag(self):
        return self._pw_280_fulldiag(Atoms('H'), hund=True, nbands=4)

    @gpwfile
    def h2_pw280_fulldiag(self):
        return self._pw_280_fulldiag(Atoms('H2', [(0, 0, 0), (0, 0, 0.7413)]),
                                     nbands=8)

    def _pw_280_fulldiag(self, atoms, **kwargs):
        atoms.set_pbc(True)
        atoms.set_cell((2., 2., 3.))
        atoms.center()
        calc = GPAW(mode=PW(280, force_complex_dtype=True),
                    xc='LDA',
                    basis='dzp',
                    parallel={'domain': 1},
                    convergence={'density': 1.e-6},
                    **kwargs)
        atoms.calc = calc
        atoms.get_potential_energy()
        calc.diagonalize_full_hamiltonian(nbands=80)
        return calc

    @gpwfile
    def fe_pw_distorted(self):
        xc = 'revTPSS'
        m = [2.9]
        fe = bulk('Fe')
        fe.set_initial_magnetic_moments(m)
        k = 3
        fe.calc = GPAW(mode=PW(800),
                       h=0.15,
                       occupations=FermiDirac(width=0.03),
                       xc=xc,
                       kpts=(k, k, k),
                       convergence={'energy': 1e-8},
                       parallel={'domain': 1, 'augment_grids': True},
                       txt=self.path / 'fe_pw_distorted.txt')
        fe.set_cell(np.dot(fe.cell,
                           [[1.02, 0, 0.03],
                            [0, 0.99, -0.02],
                            [0.2, -0.01, 1.03]]),
                    scale_atoms=True)
        fe.get_potential_energy()
        return fe.calc

    @gpwfile
    def si_pw_distorted(self):
        xc = 'TPSS'
        si = bulk('Si')
        k = 3
        si.calc = GPAW(mode=PW(250),
                       mixer=Mixer(0.7, 5, 50.0),
                       xc=xc,
                       occupations=FermiDirac(0.01),
                       kpts=(k, k, k),
                       convergence={'energy': 1e-8},
                       parallel={'domain': min(2, world.size)},
                       txt=self.path / 'si_pw_distorted.txt')
        si.set_cell(np.dot(si.cell,
                           [[1.02, 0, 0.03],
                            [0, 0.99, -0.02],
                            [0.2, -0.01, 1.03]]),
                    scale_atoms=True)
        si.get_potential_energy()
        return si.calc


@pytest.fixture(scope='session', params=sorted(_all_gpw_methodnames))
def all_gpw_files(request, gpw_files, pytestconfig):
    """This fixture parametrizes a test over all gpw_files.

    For example pytest test_generate_gpwfiles.py -n 16 is a way to quickly
    generate all gpw files independently of the rest of the test suite."""

    # Note: Parametrizing over _all_gpw_methodnames must happen *after*
    # it is populated, i.e., further down in the file than
    # the @gpwfile decorator.

    import os
    gpaw_new = os.environ.get('GPAW_NEW')

    # TODO This xfail-information should probably live closer to the
    # gpwfile definitions and not here in the fixture.
    skip_if_new = {'Cu3Au_qna', 'nicl2_pw', 'v2br4_pw_nosym', 'v2br4_pw'}
    if gpaw_new and request.param in skip_if_new:
        pytest.xfail(f'{request.param} gpwfile not yet working with GPAW_NEW')

    # Accessing each file via __getitem__ executes the calculation:
    return gpw_files[request.param]


class GPAWPlugin:
    def __init__(self):
        if world.rank == -1:
            print()
            info()

    def pytest_terminal_summary(self, terminalreporter, exitstatus, config):
        from gpaw.mpi import size
        terminalreporter.section('GPAW-MPI stuff')
        terminalreporter.write(f'size: {size}\n')


@pytest.fixture
def sg15_hydrogen():
    from io import StringIO
    from gpaw.test.pseudopotential.H_sg15 import pp_text
    from gpaw.upf import read_sg15
    # We can't easily load a non-python file from the test suite.
    # Therefore we load the pseudopotential from a Python file.
    return read_sg15(StringIO(pp_text))


def pytest_configure(config):
    # Allow for fake cupy:
    os.environ['GPAW_CPUPY'] = '1'

    if world.rank != 0:
        try:
            tw = config.get_terminal_writer()
        except AttributeError:
            pass
        else:
            tw._file = devnull
    config.pluginmanager.register(GPAWPlugin(), 'pytest_gpaw')
    for line in [
        'ci: test included in CI',
        'do: Direct optimization',
        'dscf: Delta-SCF',
        'elph: Electron-phonon',
        'fast: fast test',
        'generate_gpw_files: Dummy test to trigger gpw file precalculation',
        'gllb: GLLBSC tests',
        'gpu: GPU test',
        'hybrids: Hybrid functionals',
        'intel: fails on INTEL toolchain',
        'kspair: tests of kspair in the response code',
        'later: know failure for new refactored GPAW',
        'legacy: Old stuff that will be removed later',
        'libxc: LibXC requirered',
        'lrtddft: Linear-response TDDFT',
        'mgga: MGGA test',
        'mom: MOM',
        'ofdft: Orbital-free DFT',
        'response: tests of the response code',
        'rpa: tests of RPA',
        'rttddft: Real-time TDDFT',
        'serial: run in serial only',
        'sic: PZ-SIC',
        'slow: slow test',
        'soc: Spin-orbit coupling',
        'stress: Calculation of stress tensor',
        'wannier: Wannier functions',
        'pipekmezey : PipekMezey wannier functions']:
        config.addinivalue_line('markers', line)


def pytest_runtest_setup(item):
    """Skip some tests.

    If:

    * they depend on libxc and GPAW is not compiled with libxc
    * they are before $PYTEST_START_AFTER
    """
    from gpaw import libraries

    if world.size > 1:
        for mark in item.iter_markers():
            if mark.name == 'serial':
                pytest.skip('Only run in serial')

    if item.location[0] <= os.environ.get('PYTEST_START_AFTER', ''):
        pytest.skip('Not after $PYTEST_START_AFTER')
        return

    if libraries['libxc']:
        return

    if any(mark.name in {'libxc', 'mgga'}
           for mark in item.iter_markers()):
        pytest.skip('No LibXC.')


@pytest.fixture
def scalapack():
    """Skip if not compiled with sl.

    This fixture otherwise does not return or do anything."""
    from gpaw.utilities import compiled_with_sl
    if not compiled_with_sl():
        pytest.skip('no scalapack')


@pytest.fixture
def needs_ase_master():
    from ase.utils.filecache import MultiFileJSONCache
    try:
        MultiFileJSONCache('bla-bla', comm=None)
    except TypeError:
        pytest.skip('ASE is too old')


def pytest_report_header(config, startdir):
    # Use this to add custom information to the pytest printout.
    yield f'GPAW MPI rank={world.rank}, size={world.size}'

    # We want the user to be able to see where gpw files are cached,
    # but the only way to see the cache location is to make a directory
    # inside it.  mkdir('') returns the toplevel cache dir without
    # actually creating a subdirectory:
    cachedir = config.cache.mkdir('')
    yield f'Cache directory including gpw files: {cachedir}'


@pytest.fixture
def rng():
    """Seeded random number generator.

    Tests should be deterministic and should use this
    fixture or initialize their own rng."""
    return np.random.default_rng(42)


@pytest.fixture
def gpaw_new() -> bool:
    """Are we testing the new code?"""
    return os.environ.get('GPAW_NEW')
