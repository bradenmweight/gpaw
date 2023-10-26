import pytest
import numpy as np

from gpaw.core import UGDesc
from gpaw.gpu import cupy as cp
from gpaw.mpi import world
from gpaw.spline import Spline


@pytest.mark.parametrize('dtype', [float, complex])
def test_lfc(dtype):
    s = Spline(0, 1.0, [1.0, 0.5, 0.0])
    n = 40
    a = 8.0
    fracpos_ac = [(0.5, 0.5, 0.25 + 0.25 * i) for i in [0, 1, 2]]

    grid = UGDesc(cell=[a, a, a], size=(n, n, n), comm=world, dtype=dtype)
    basis_cpu = grid.atom_centered_functions([[s],[s],[s]],
                                             positions=fracpos_ac, xp=np)
    basis_gpu = grid.atom_centered_functions([[s],[s],[s]],
                                             positions=fracpos_ac, xp=cp)

    P_cpu_ani = basis_cpu.layout.empty()
    P_gpu_ani = basis_gpu.layout.empty()
    P_cpu_ani.data[:] = 1.0
    P_gpu_ani.data[:] = 1.0

    b_cpu = grid.zeros(xp=np)
    b_gpu = grid.zeros(xp=cp)
    basis_cpu.add_to(b_cpu, P_cpu_ani)
    basis_gpu.add_to(b_gpu, P_gpu_ani)

    assert b_cpu.data == pytest.approx(b_gpu.data.get(), abs=1e-12)

    out_cpu_ani = basis_cpu.integrate(b_cpu)
    out_gpu_ani = basis_gpu.integrate(b_gpu)

    for a, out_cpu_ni in out_cpu_ani.items():
        assert out_cpu_ni == pytest.approx(out_gpu_ani[a].get(), abs=1e-12)
