import numpy as np
from typing import TypeVar, Any
def adjust_positions() -> None: ...
def adjust_momenta() -> None: ...
def calculate_forces_H2O() -> None: ...
def localize(Z_nnc: Any, U_nn: Any) -> float: ...
def spherical_harmonics() -> None: ...
def hartree(l: int, nrdr: np.ndarray, r: np.ndarray, vr: np.ndarray) -> None: ...
def get_num_threads() -> int: ...
def pack(A: np.ndarray) -> np.ndarray: ...
T = TypeVar('T', float, complex)
def mmm(alpha: T,
        a: np.ndarray,
        opa: str,
        b: np.ndarray,
        opb: str,
        beta: T,
        c: np.ndarray) -> None: ...
def gemm(alpha, a, b, beta, c, transa='n') -> None: ...
def rk(alpha, a, beta, c, trans='c') -> None: ...
def r2k(alpha, a, b, beta, c, trans='c') -> None: ...
def gemmdot() -> None: ...
class Communicator:
    rank: int
    size: int
def add_to_density(f: float, psit: np.ndarray, density: np.ndarray) -> None: ...
