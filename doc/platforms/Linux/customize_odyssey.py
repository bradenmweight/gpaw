scalapack = True

mklpath ='/n/sw/intel/mkl/10.3.1.107/composerxe-2011.3.174/mkl'
omppath ='/n/sw/openmpi-1.5.3_intel-12.3.174ib'

compiler = 'icc'

libraries = ['mpi', 'mpi_f77', 'mkl_scalapack_lp64', 'mkl_lapack95_lp64', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_mc', 'mkl_core', 'mkl_def', 'mkl_intel_thread', 'iomp5']
library_dirs += [f'{omppath}/lib', f'{mklpath}/lib/intel64']
include_dirs += ['/usr/include', f'{omppath}/include', f'{mklpath}/include']

extra_link_args += [f'{mklpath}/lib/intel64/libmkl_blacs_openmpi_lp64.a', f'{mklpath}/lib/intel64/libmkl_blas95_lp64.a']

extra_compile_args += ['-O3', '-std=c99', '-w']

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

mpicompiler = 'mpicc'
mpilinker = mpicompiler
