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
def pw_insert(coef_G: np.ndarray, Q_G: np.ndarray, s: float, array_Q: np.ndarray) -> None: ...
def pblas_tran(N: int, M: int,
               alpha: float, a_MN: np.ndarray,
               beta:float, c_NM: np.ndarray,
               desca: np.ndarray, descc: np.ndarray,
               conj: bool) -> None: ...
def scalapack_set(a: np.ndarray, desc: np.ndarray, alpha: float, beta: float,
                  uplo: str, n: int, m: int, ja: int, ia: int) -> None: ...
