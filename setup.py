#!/usr/bin/env python
# Copyright (C) 2003-2020  CAMP
# Please see the accompanying LICENSE file for further information.

import os
import re
import sys
import warnings
from pathlib import Path
from subprocess import run
from sysconfig import get_platform

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install

from config import build_interpreter, check_dependencies, write_configuration

assert sys.version_info >= (3, 7)


def warn_deprecated(msg, error=False):
    msg = f'\n\n{msg}\n\n'
    if error:
        raise ValueError(msg)
    else:
        warnings.warn(msg, DeprecationWarning)


# Get the current version number:
txt = Path('gpaw/__init__.py').read_text()
version = re.search("__version__ = '(.*)'", txt)[1]
ase_version_required = re.search("__ase_version_required__ = '(.*)'", txt)[1]

description = 'GPAW: DFT and beyond within the projector-augmented wave method'
long_description = Path('README.rst').read_text()

for i, arg in enumerate(sys.argv):
    if arg.startswith('--customize='):
        custom = arg.split('=')[1]
        raise DeprecationWarning(
            f'Please set GPAW_CONFIG={custom} or place {custom} in ' +
            '~/.gpaw/siteconfig.py')

libraries = ['xc']
library_dirs = []
include_dirs = []
extra_link_args = []
extra_compile_args = ['-Wall', '-Wno-unknown-pragmas', '-std=c99']
runtime_library_dirs = []
extra_objects = []
define_macros = [('NPY_NO_DEPRECATED_API', '7'),
                 ('GPAW_NO_UNDERSCORE_CBLACS', '1'),
                 ('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
if os.getenv('GPAW_GPU'):
    define_macros.append(('GPAW_GPU_AWARE_MPI', '1'))
undef_macros = ['NDEBUG']

parallel_python_interpreter = False
compiler = None
mpi = False
noblas = False
nolibxc = False
fftw = False
scalapack = False
libvdwxc = False
elpa = False

# Advanced:
# If these are defined, they replace
# all the default args from setuptools
compiler_args = None
linker_so_args = None
linker_exe_args = None

# MPI is enabled by default if `mpicc` is found
found_mpicc = (os.name != 'nt'
               and run(['which', 'mpicc'],
                       capture_output=True).returncode == 0)
if found_mpicc:
    mpi = True

# Search and store current git hash if possible
try:
    from ase.utils import search_current_git_hash
    githash = search_current_git_hash('gpaw')
    if githash is not None:
        define_macros += [('GPAW_GITHASH', githash)]
    else:
        print('.git directory not found. GPAW git hash not written.')
except ImportError:
    print('ASE not found. GPAW git hash not written.')

# User provided customizations:
gpaw_config = os.environ.get('GPAW_CONFIG')
if gpaw_config and not Path(gpaw_config).is_file():
    raise FileNotFoundError(gpaw_config)
for siteconfig in [gpaw_config,
                   'siteconfig.py',
                   '~/.gpaw/siteconfig.py']:
    if siteconfig is not None:
        path = Path(siteconfig).expanduser()
        if path.is_file():
            print('Reading configuration from', path)
            exec(path.read_text())
            break
else:  # no break
    if not noblas:
        libraries.append('blas')

if 'mpicompiler' in locals():
    mpicompiler = locals()['mpicompiler']
    msg = 'Please remove deprecated declaration of mpicompiler.'
    if mpicompiler is None:
        mpi = False
        msg += (' Define instead in siteconfig one of the following lines:'
                '\n\nmpi = False\nmpi = True')
    else:
        mpi = True
        compiler = mpicompiler
        msg += (' Define instead in siteconfig:'
                f'\n\nmpi = True\ncompiler = {repr(compiler)}')
    warn_deprecated(msg)

if 'mpilinker' in locals():
    mpilinker = locals()['mpilinker']
    msg = ('Please remove deprecated declaration of mpilinker:'
           f'\ncompiler={repr(compiler)} will be used for linking.')
    if mpilinker == compiler:
        warn_deprecated(msg)
    else:
        msg += ('\nPlease contact GPAW developers if you need '
                'different commands for linking and compiling.')
        warn_deprecated(msg, error=True)

for key in ['libraries', 'library_dirs', 'include_dirs',
            'runtime_library_dirs', 'define_macros']:
    mpi_key = 'mpi_' + key
    if mpi_key in locals():
        warn_deprecated(
            f'Please remove deprecated declaration of {mpi_key}'
            ' and use only {key} instead.'
            f'\nAdding {mpi_key} to {key}.')
        locals()[key] += locals()[mpi_key]

if mpi:
    if compiler is None:
        if found_mpicc:
            compiler = 'mpicc'
        else:
            raise ValueError('Define compiler for MPI in siteconfig:'
                             "\n\ncompiler = '...'")

if parallel_python_interpreter:
    parallel_python_exefile = None
    if not mpi:
        raise ValueError('MPI is needed for parallel_python_interpreter.'
                         ' Define in siteconfig:'
                         '\nparallel_python_interpreter = True'
                         '\nmpi = True'
                         "\ncompiler = '...'  # MPI compiler, e.g., 'mpicc'"
                         )

platform_id = os.getenv('CPU_ARCH')
if platform_id:
    os.environ['_PYTHON_HOST_PLATFORM'] = get_platform() + '-' + platform_id

for flag, name in [(noblas, 'GPAW_WITHOUT_BLAS'),
                   (nolibxc, 'GPAW_WITHOUT_LIBXC'),
                   (mpi, 'PARALLEL'),
                   (fftw, 'GPAW_WITH_FFTW'),
                   (scalapack, 'GPAW_WITH_SL'),
                   (libvdwxc, 'GPAW_WITH_LIBVDWXC'),
                   (elpa, 'GPAW_WITH_ELPA')]:
    if flag:
        define_macros.append((name, '1'))

sources = [Path('c/bmgs/bmgs.c')]
sources += Path('c').glob('*.c')
sources += Path('c/xc').glob('*.c')
if nolibxc:
    for name in ['libxc.c', 'm06l.c',
                 'tpss.c', 'revtpss.c', 'revtpss_c_pbe.c',
                 'xc_mgga.c']:
        sources.remove(Path(f'c/xc/{name}'))
    if 'xc' in libraries:
        libraries.remove('xc')

# Make build process deterministic (for "reproducible build")
sources = [str(source) for source in sources]
sources.sort()

check_dependencies(sources)

# Convert Path objects to str:
library_dirs = [str(dir) for dir in library_dirs]
include_dirs = [str(dir) for dir in include_dirs]

extensions = [Extension('_gpaw',
                        sources,
                        libraries=libraries,
                        library_dirs=library_dirs,
                        include_dirs=include_dirs,
                        define_macros=define_macros,
                        undef_macros=undef_macros,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args,
                        runtime_library_dirs=runtime_library_dirs,
                        extra_objects=extra_objects,
                        language='c')]


if os.environ.get('GPAW_GPU'):
    # Hardcoded for LUMI right now!
    target = os.environ.get('HCC_AMDGPU_TARGET', 'gfx90a')
    # TODO: Build this also via extension
    assert os.system(
        f'HCC_AMDGPU_TARGET={target} hipcc -fPIC -fgpu-rdc '
        '-c c/gpu/hip_kernels.cpp -o c/gpu/hip_kernels.o') == 0
    assert os.system(
        f'HCC_AMDGPU_TARGET={target} hipcc -shared -fgpu-rdc --hip-link '
        '-o c/gpu/hip_kernels.so c/gpu/hip_kernels.o') == 0

    extensions.append(
        Extension('_gpaw_gpu',
                  ['c/gpu/gpaw_gpu.c'],
                  libraries=[],
                  library_dirs=['c/gpu'],
                  setup_requires=['numpy'],
                  include_dirs=include_dirs,
                  define_macros=[('NPY_NO_DEPRECATED_API', 7)],
                  undef_macros=[],
                  extra_link_args=[
                      f'-Wl,-rpath={Path("c/gpu").resolve()}'],
                  extra_compile_args=['-std=c99'],
                  # ,'-Werror=implicit-function-declaration'],
                  runtime_library_dirs=['c/gpu'],
                  extra_objects=[
                      str(Path('c/gpu/hip_kernels.so').resolve())]))


write_configuration(define_macros, include_dirs, libraries, library_dirs,
                    extra_link_args, extra_compile_args,
                    runtime_library_dirs, extra_objects, compiler)


class build_ext(_build_ext):
    def run(self):
        import numpy as np
        self.include_dirs.append(np.get_include())
        super().run()

    def build_extensions(self):
        # Override the compiler executables
        # A hack to change the used compiler and linker, inspired by
        # https://shwina.github.io/custom-compiler-linker-extensions/
        for (name, my_args) in [('compiler', compiler_args),
                                ('compiler_so', compiler_args),
                                ('linker_so', linker_so_args),
                                ('linker_exe', linker_exe_args)]:
            new_args = []
            old_args = getattr(self.compiler, name)
            # Set executable
            if compiler is not None:
                new_args += [compiler]
            else:
                new_args += [old_args[0]]
            # Set args
            if my_args is not None:
                new_args += my_args
            else:
                new_args += old_args[1:]
            self.compiler.set_executable(name, new_args)

        super().build_extensions()

        if parallel_python_interpreter:
            global parallel_python_exefile

            assert len(self.extensions) == 1, \
                'Fix gpaw-python build for multiple extensions'
            extension = self.extensions[0]

            # Path for the bin (analogous to build_lib)
            build_bin = Path(str(self.build_lib).replace('lib', 'bin'))

            # List of object files already built for the extension
            objects = []
            for src in sources:
                obj = Path(self.build_temp) / Path(src).with_suffix('.o')
                objects.append(str(obj))

            # Build gpaw-python
            parallel_python_exefile = build_interpreter(
                self.compiler, extension, objects,
                build_temp=self.build_temp,
                build_bin=build_bin,
                debug=self.debug)

        print("Build temp:", self.build_temp)
        print("Build lib: ", self.build_lib)
        if parallel_python_interpreter:
            print("Build bin: ", build_bin)


class install(_install):
    def run(self):
        super().run()
        if parallel_python_interpreter:
            self.copy_file(parallel_python_exefile, self.install_scripts)


class develop(_develop):
    def run(self):
        super().run()
        if parallel_python_interpreter:
            self.copy_file(parallel_python_exefile, self.script_dir)


cmdclass = {'build_ext': build_ext,
            'install': install,
            'develop': develop}

files = ['gpaw-analyse-basis', 'gpaw-basis',
         'gpaw-plot-parallel-timings', 'gpaw-runscript',
         'gpaw-setup', 'gpaw-upfplot']
scripts = [str(Path('tools') / script) for script in files]


setup(name='gpaw',
      version=version,
      description=description,
      long_description=long_description,
      maintainer='GPAW-community',
      maintainer_email='gpaw-users@listserv.fysik.dtu.dk',
      url='https://wiki.fysik.dtu.dk/gpaw',
      license='GPLv3+',
      platforms=['unix'],
      packages=find_packages(),
      package_data={'gpaw': ['py.typed']},
      entry_points={
          'console_scripts': ['gpaw = gpaw.cli.main:main']},
      setup_requires=['numpy'],
      install_requires=[f'ase>={ase_version_required}',
                        'scipy>=1.2.0',
                        'pyyaml'],
      extras_require={'docs': ['sphinx-rtd-theme',
                               'graphviz'],
                      'devel': ['flake8',
                                'mypy',
                                'pytest-xdist',
                                'interrogate']},
      ext_modules=extensions,
      scripts=scripts,
      cmdclass=cmdclass,
      classifiers=[
          'Development Status :: 6 - Mature',
          'License :: OSI Approved :: '
          'GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Topic :: Scientific/Engineering :: Physics'])
