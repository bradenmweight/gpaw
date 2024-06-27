#!/usr/bin/env python3
"""Install gpaw on Niflheim in a virtual environment.

Also installs ase, ase-ext, spglib, sklearn and myqueue.
"""
import argparse
import os
import shutil
import subprocess
from pathlib import Path
from sys import version_info

if version_info < (3, 8):
    raise ValueError('Please use Python-3.8 or later')

# Python version in the venv that we are creating
version = '3.11'
fversion = 'cpython-311'

# Niflheim login hosts, with the oldest architecture as the first
nifllogin = ['slid2',  # broadwell_el8 (xeon24el8)
             'thul',  # skylake_el8 (xeon40el8)
             'surt',  # icelake (xeon56)
             'fjorm']  # epyc9004 (epyc96)


# Easybuild uses a hierarchy of toolchains for the main foss and intel
# chains.  The order in the tuples before are
#  fullchain: Full chain.
#  mathchain: Chain with math libraries but no MPI
#  compchain: Chain with full compiler suite (but no fancy libs)
#  corechain: Core compiler
# The subchain complementary to 'mathchain', with MPI but no math libs, is
# not used here.

_gcccore = 'GCCcore-12.3.0'
toolchains = {
    'foss': dict(
        fullchain='foss-2023a',
        mathchain='gfbf-2023a',
        compchain='GCC-12.3.0',
        corechain=_gcccore,
    ),
    'intel': dict(
        fullchain='intel-2023a',
        mathchain='iimkl-2023a',
        compchain='intel-compilers-2023.1.0',
        corechain=_gcccore,
    )
}

# These modules are always loaded
module_cmds_all = """\
module purge
unset PYTHONPATH
module load GPAW-setups/24.1.0
module load ELPA/2023.05.001-{fullchain}
module load Wannier90/3.1.0-{fullchain}
module load Python-bundle-PyPI/2023.06-{corechain}
module load Tkinter/3.11.3-{corechain}
module load libxc/6.2.2-{compchain}
"""

# These modules are not loaded if --piponly is specified
module_cmds_easybuild = """\
module load matplotlib/3.7.2-{mathchain}
module load scikit-learn/1.3.1-{mathchain}
module load spglib-python/2.1.0-{mathchain}
"""

# These modules are loaded depending on the toolchain
module_cmds_tc = {
    'foss': """\
module load libvdwxc/0.4.0-{fullchain}
""",
    'intel': ""
}

module_cmds_arch_dependent = """\
if [ "$CPU_ARCH" == "icelake" ];\
then module load CuPy/12.3.0-{fullchain}-CUDA-12.1.1;fi
"""


activate_extra = """
export GPAW_SETUP_PATH=$GPAW_SETUP_PATH:{venv}/gpaw-basis-pvalence-0.9.20000

# Set matplotlib backend:
if [[ $SLURM_SUBMIT_DIR ]]; then
    export MPLBACKEND=Agg
    export PYTHONWARNINGS="ignore:Matplotlib is currently using agg"
else
    export MPLBACKEND=TkAgg
fi
"""

dftd3 = """\
mkdir {venv}/DFTD3
cd {venv}/DFTD3
URL=https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3
wget $URL/dftd3.tgz
tar -xf dftd3.tgz
ssh {nifllogin[0]} ". {venv}/bin/activate && cd {venv}/DFTD3 && make >& d3.log"
ln -s {venv}/DFTD3/dftd3 {venv}/bin
"""


def run(cmd: str, **kwargs) -> subprocess.CompletedProcess:
    print(cmd)
    return subprocess.run(cmd, shell=True, check=True, **kwargs)


def compile_gpaw_c_code(gpaw: Path, activate: Path) -> None:
    """Compile for all architectures: xeon24, xeon40, ..."""
    # Remove targets:
    for path in gpaw.glob('build/lib.linux-x86_64-*/_gpaw.*.so'):
        path.unlink()

    # Compile:
    for host in nifllogin:
        run(f'ssh {host} ". {activate} && pip install -q -e {gpaw}"')

    # Clean up:
    for path in gpaw.glob('_gpaw.*.so'):
        path.unlink()
    for path in gpaw.glob('build/temp.linux-x86_64-*'):
        shutil.rmtree(path)


def fix_installed_scripts(venvdir: Path,
                          rootdir: str,
                          pythonroot: str) -> None:
    """Fix command line tools so they work in the virtual environment.

    Command line tools (pytest, sphinx-build etc) fail in virtual
    enviroments created with --system-site-packages, as the scripts
    are not copied into the virtual environment.  The scripts have
    the original Python interpreter hardcoded in the hash-bang line.

    This function copies all scripts into the virtual environment,
    and changes the hash-bang so it works.  Starting with the 2023a
    toolchains, the scripts are distributed over more than one
    EasyBuild module.

    Arguments:
    venvdir: Path to the virtual environment
    rootdir: string holding folder of the EasyBuild package being processed
    pythondir: string holding folder of the Python package.
    """

    assert rootdir is not None
    assert pythonroot is not None
    bindir = rootdir / Path('bin')
    print(f'Patching executable scripts from {bindir} to {venvdir}/bin')
    assert '+' not in str(pythonroot) and '+' not in str(venvdir), (
        'Script will fail with "+" in folder names!')
    sedscript = f's+{pythonroot}+{venvdir}+g'

    # Loop over potential executables
    for exe in bindir.iterdir():
        target = venvdir / 'bin' / exe.name
        # Skip files that already exist, are part of Python itself,
        # or are not a regular file or symlink to a file.
        if (not target.exists()
                and not exe.name.lower().startswith('python')
                and exe.is_file()):
            # Check if it is a script file referring the original
            # Python executable in the hash-bang
            with open(exe) as f:
                firstline = f.readline()
            if pythonroot in firstline:
                shutil.copy2(exe, target, follow_symlinks=False)
                # Now patch the file (if not a symlink)
                if not exe.is_symlink():
                    assert not target.is_symlink()
                    subprocess.run(
                        f"sed -e '{sedscript}' --in-place '{target}'",
                        shell=True,
                        check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('venv', help='Name of venv.')
    parser.add_argument('--toolchain', default='foss',
                        choices=['foss', 'intel'],
                        help='Default is foss.')
    parser.add_argument('--dftd3', action='store_true',
                        help='Also build DFTD3.')
    parser.add_argument('--recompile', action='store_true',
                        help='Recompile the GPAW C-extensions in an '
                        'exising venv.')
    parser.add_argument('--piponly', action='store_true',
                        help='Do not use EasyBuild python modules, '
                        'install from pip (may affect performance).')
    args = parser.parse_args()

    # if args.toolchain == 'intel':
    #     raise ValueError('See: https://gitlab.com/gpaw/gpaw/-/issues/241')

    venv = Path(args.venv).absolute()
    activate = venv / 'bin/activate'
    gpaw = venv / 'gpaw'

    if args.recompile:
        compile_gpaw_c_code(gpaw, activate)
        return 0

    # Sanity checks
    if args.toolchain not in ('foss', 'intel'):
        raise ValueError(f'Unsupported toolchain "{args.toolchain}"')

    module_cmds = module_cmds_all.format(**toolchains[args.toolchain])
    if not args.piponly:
        module_cmds += module_cmds_easybuild.format(
            **toolchains[args.toolchain])
    module_cmds += module_cmds_tc[args.toolchain].format(
        **toolchains[args.toolchain])
    module_cmds += module_cmds_arch_dependent.format(
        **toolchains[args.toolchain])
    cmds = (' && '.join(module_cmds.splitlines()) +
            f' && python3 -m venv --system-site-packages {args.venv}')
    run(cmds)

    os.chdir(venv)

    activate.write_text(module_cmds +
                        activate.read_text())

    run(f'. {activate} && pip install --upgrade pip -q')

    # Fix venv so pytest etc work
    pythonroot = None
    for ebrootvar in ('EBROOTPYTHON', 'EBROOTPYTHONMINBUNDLEMINPYPI'):
        # Note that we need the environment variable from the newly
        # created venv, NOT from this process!
        comm = run(f'. {activate} && echo ${ebrootvar}',
                   capture_output=True, text=True)
        ebrootdir = comm.stdout.strip()
        if pythonroot is None:
            # The first module is the actual Python module.
            pythonroot = ebrootdir
        assert ebrootdir, f'Env variable {ebrootvar} appears to be unset.'
        fix_installed_scripts(venvdir=venv,
                              rootdir=ebrootdir,
                              pythonroot=pythonroot)

    packages = ['myqueue',
                'graphviz',
                'qeh',
                'sphinx_rtd_theme']
    if args.piponly:
        packages += ['matplotlib',
                     'scipy',
                     'pandas',
                     'pytest',
                     'pytest-xdist',
                     'pytest-mock',
                     'scikit-learn']
    run(f'. {activate} && pip install -q -U ' + ' '.join(packages))

    for name in ['ase', 'gpaw']:
        run(f'git clone -q https://gitlab.com/{name}/{name}.git')

    run(f'. {activate} && pip install -q -e ase/')

    if args.dftd3:
        run(' && '.join(dftd3.format(venv=venv,
                                     nifllogin=nifllogin).splitlines()))

    # Compile ase-ext C-extension on old thul so that it works on
    # newer architectures
    run(f'ssh {nifllogin[0]} ". {activate} && pip install -q ase-ext"')

    if args.piponly:
        run('git clone -q https://github.com/spglib/spglib.git')
        run(f'ssh {nifllogin[0]} ". {activate} && pip install {venv}/spglib"')

    # Install GPAW:
    siteconfig = Path(
        f'gpaw/doc/platforms/Linux/Niflheim/siteconfig-{args.toolchain}.py')
    Path('gpaw/siteconfig.py').write_text(siteconfig.read_text())

    compile_gpaw_c_code(gpaw, activate)

    for fro, to in [('ivybridge', 'sandybridge'),
                    ('nahelem', 'icelake')]:
        f = gpaw / f'build/lib.linux-x86_64-{fro}-{fversion}'
        t = gpaw / f'build/lib.linux-x86_64-{to}-{fversion}'
        f.symlink_to(t)

    # Create .pth file to load correct .so file:
    pth = (
        'import sys, os; '
        'arch = os.environ["CPU_ARCH"]; '
        f"path = f'{venv}/gpaw/build/lib.linux-x86_64-{{arch}}-{fversion}'; "
        'sys.path.append(path)\n')
    Path(f'lib/python{version}/site-packages/niflheim.pth').write_text(pth)

    # Install extra basis-functions:
    run(f'. {activate} && gpaw install-data --basis --version=20000 '
        f'{venv} --no-register')

    extra = activate_extra.format(venv=venv)

    # Tab completion:
    for cmd in ['ase', 'gpaw', 'mq', 'pip']:
        if cmd == 'gpaw':
            # Currently, running the "gpaw" command writes warning message
            # to stdout, so "gpaw completion" does not work!
            continue
        txt = run(f'. {activate} && {cmd} completion' +
                  (' --bash' if cmd == 'pip' else ''),
                  capture_output=True).stdout.decode()
        extra += txt

    # gpaw-hack:
    python = venv / 'bin/python3'
    complete = venv / 'gpaw/gpaw/cli/complete.py'
    extra += f'complete -o default -C "{python} {complete}" gpaw\n'

    activate.write_text(activate.read_text() + extra)

    # Run tests:
    run(f'. {activate} && ase info && gpaw test')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
