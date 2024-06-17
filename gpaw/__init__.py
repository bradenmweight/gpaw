# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.
"""Main gpaw module."""
from __future__ import annotations
import os
import sys
import contextlib
from pathlib import Path
from typing import List, Union, Any, TYPE_CHECKING

__version__ = '24.7.0b1'
__ase_version_required__ = '3.23.0'

__all__ = ['GPAW',
           'Mixer', 'MixerSum', 'MixerDif', 'MixerSum2',
           'MixerFull',
           'CG', 'Davidson', 'RMMDIIS', 'DirectLCAO',
           'PoissonSolver',
           'FermiDirac', 'MethfesselPaxton', 'MarzariVanderbilt',
           'PW', 'LCAO', 'FD',
           'restart']

setup_paths: List[Union[str, Path]] = []
is_gpaw_python = '_gpaw' in sys.builtin_module_names
dry_run = 0

# When type-checking or running pytest, we want the debug-wrappers enabled:
debug: bool = (TYPE_CHECKING or
               'pytest' in sys.modules or
               bool(sys.flags.debug))


@contextlib.contextmanager
def disable_dry_run():
    """Context manager for temporarily disabling dry-run mode.

    Useful for skipping exit in the GPAW constructor.
    """
    global dry_run
    size = dry_run
    dry_run = 0
    yield
    dry_run = size


def get_scipy_version():
    import scipy
    # This is in a function because we don't like to have the scipy
    # import at module level
    return [int(x) for x in scipy.__version__.split('.')[:2]]


if 'OMP_NUM_THREADS' not in os.environ:
    os.environ['OMP_NUM_THREADS'] = '1'


class ConvergenceError(Exception):
    pass


class KohnShamConvergenceError(ConvergenceError):
    pass


class PoissonConvergenceError(ConvergenceError):
    pass


class KPointError(Exception):
    pass


class BadParallelization(Exception):
    """Error indicating missing parallelization support."""
    pass


def get_libraries() -> dict[str, str]:
    import gpaw.cgpaw as cgpaw
    libraries: dict[str, str] = {}
    if hasattr(cgpaw, 'lxcXCFunctional'):
        libraries['libxc'] = getattr(cgpaw, 'libxc_version', '2.x.y')
    else:
        libraries['libxc'] = ''
    return libraries


def parse_arguments(argv):
    from argparse import (ArgumentParser, REMAINDER,
                          RawDescriptionHelpFormatter)
    import warnings

    # With gpaw-python BLAS symbols are in global scope and we need to
    # ensure that NumPy and SciPy use symbols from their own dependencies
    if is_gpaw_python:
        old_dlopen_flags = sys.getdlopenflags()
        sys.setdlopenflags(old_dlopen_flags | os.RTLD_DEEPBIND)

    if is_gpaw_python:
        sys.setdlopenflags(old_dlopen_flags)

    version = sys.version.replace('\n', '')
    p = ArgumentParser(usage='%(prog)s [OPTION ...] [-c | -m] SCRIPT'
                       ' [ARG ...]',
                       description='Run a parallel GPAW calculation.\n\n'
                       f'Compiled with:\n  Python {version}',
                       formatter_class=RawDescriptionHelpFormatter)

    p.add_argument('--command', '-c', action='store_true',
                   help='execute Python string given as SCRIPT')
    p.add_argument('--module', '-m', action='store_true',
                   help='run library module given as SCRIPT')
    p.add_argument('-W', metavar='argument',
                   action='append', default=[], dest='warnings',
                   help='warning control.  See the documentation of -W for '
                   'the Python interpreter')
    p.add_argument('script', metavar='SCRIPT',
                   help='calculation script')
    p.add_argument('options', metavar='ARG',
                   help='arguments forwarded to SCRIPT', nargs=REMAINDER)

    args = p.parse_args(argv[1:])

    if args.command and args.module:
        p.error('-c and -m are mutually exclusive')

    sys.argv = [args.script] + args.options

    for w in args.warnings:
        # Need to convert between python -W syntax to call
        # warnings.filterwarnings():
        warn_args = w.split(':')
        assert len(warn_args) <= 5

        if warn_args[0] == 'all':
            warn_args[0] = 'always'
        if len(warn_args) >= 3:
            # e.g. 'UserWarning' (string) -> UserWarning (class)
            warn_args[2] = globals().get(warn_args[2])
        if len(warn_args) == 5:
            warn_args[4] = int(warn_args[4])

        warnings.filterwarnings(*warn_args, append=True)

    return args


def __getattr__(attr: str) -> Any:
    """
    Implement the lazy importing of classes in submodules."""

    import importlib

    try:
        import_target = all_lazy_imports[attr]
    except KeyError:
        raise AttributeError(
            f'{__getattr__.__module__}: no attribute named `.{attr!r}`'
        ) from None

    module, sep, target = import_target.rpartition('.')
    assert module and all(chunk.isidentifier() for chunk in module.split('.'))
    assert sep
    assert target.isidentifier()
    return getattr(importlib.import_module(module), target)


def __dir__() -> List[str]:
    """
    Get the normally-present module attributes and the lazily-imported objects.
    """
    return [*globals(), *all_lazy_imports]


all_lazy_imports = dict(
    OldGPAW='gpaw.calculator.GPAW',
    Mixer='gpaw.mixer.Mixer',
    MixerSum='gpaw.mixer.MixerSum',
    MixerDif='gpaw.mixer.MixerDif',
    MixerSum2='gpaw.mixer.MixerSum2',
    MixerFull='gpaw.mixer.MixerFull',

    Davidson='gpaw.eigensolvers.Davidson',
    RMMDIIS='gpaw.eigensolvers.RMMDIIS',
    CG='gpaw.eigensolvers.CG',
    DirectLCAO='gpaw.eigensolvers.DirectLCAO',

    PoissonSolver='gpaw.poisson.PoissonSolver',
    FermiDirac='gpaw.occupations.FermiDirac',
    MethfesselPaxton='gpaw.occupations.MethfesselPaxton',
    MarzariVanderbilt='gpaw.occupations.MarzariVanderbilt',
    FD='gpaw.wavefunctions.fd.FD',
    LCAO='gpaw.wavefunctions.lcao.LCAO',
    PW='gpaw.wavefunctions.pw.PW',
)


class BroadcastImports:
    def __enter__(self):
        from gpaw._broadcast_imports import broadcast_imports
        self._context = broadcast_imports
        return self._context.__enter__()

    def __exit__(self, *args):
        self._context.__exit__(*args)


broadcast_imports = BroadcastImports()


def main():
    with broadcast_imports:
        import runpy

        # Apparently we need the scipy.linalg import for compatibility?
        import scipy.linalg  # noqa: F401

        for attr in all_lazy_imports:
            __getattr__(attr)

        gpaw_args = parse_arguments(sys.argv)
    # The normal Python interpreter puts . in sys.path, so we also do that:
    sys.path.insert(0, '.')
    # Stacktraces can be shortened by running script with
    # PyExec_AnyFile and friends.  Might be nicer
    if gpaw_args.command:
        d = {'__name__': '__main__'}
        exec(gpaw_args.script, d, d)
    elif gpaw_args.module:
        # Python has: python [-m MOD] [-c CMD] [SCRIPT]
        # We use a much better way: gpaw-python [-m | -c] SCRIPT
        runpy.run_module(gpaw_args.script, run_name='__main__')
    else:
        runpy.run_path(gpaw_args.script, run_name='__main__')


if debug:
    import numpy as np
    np.seterr(over='raise', divide='raise', invalid='raise', under='ignore')
    oldempty = np.empty
    oldempty_like = np.empty_like

    def empty(*args, **kwargs):
        a = oldempty(*args, **kwargs)
        try:
            a.fill(np.nan)
        except ValueError:
            a.fill(42)
        return a

    def empty_like(*args, **kwargs):
        a = oldempty_like(*args, **kwargs)
        try:
            a.fill(np.nan)
        except ValueError:
            a.fill(-42)
        return a

    np.empty = empty
    np.empty_like = empty_like


def GPAW(*args, **kwargs) -> Any:
    if int(os.environ.get('GPAW_NEW', '0')):
        from gpaw.new.ase_interface import GPAW as NewGPAW
        return NewGPAW(*args, **kwargs)

    from gpaw.calculator import GPAW as OldGPAW
    return OldGPAW(*args, **kwargs)


def restart(filename, Class=None, **kwargs):
    if Class is None:
        from gpaw import GPAW as Class
    calc = Class(filename, **kwargs)
    atoms = calc.get_atoms()
    return atoms, calc


def read_rc_file():
    home = os.environ.get('HOME')
    if home is not None:
        rc = os.path.join(home, '.gpaw', 'rc.py')
        if os.path.isfile(rc):
            # Read file in ~/.gpaw/rc.py
            with open(rc) as fd:
                exec(fd.read())


def initialize_data_paths():
    try:
        setup_paths[:0] = os.environ['GPAW_SETUP_PATH'].split(os.pathsep)
    except KeyError:
        if len(setup_paths) == 0:
            if os.pathsep == ';':
                setup_paths[:] = [r'C:\gpaw-setups']
            else:
                setup_paths[:] = ['/usr/local/share/gpaw-setups',
                                  '/usr/share/gpaw-setups']


read_rc_file()
initialize_data_paths()


def RMM_DIIS(*args, **kwargs):
    import warnings
    from gpaw import RMMDIIS
    warnings.warn('Please use RMMDIIS instead of RMM_DIIS')
    return RMMDIIS(*args, **kwargs)
