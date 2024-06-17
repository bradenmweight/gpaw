import numpy as np
from functools import cached_property


class ResponseKPointGrid:
    def __init__(self, kd):
        self.kd = kd

    @cached_property
    def kptfinder(self):
        from gpaw.response.symmetry import KPointFinder
        return KPointFinder(self.kd.bzk_kc)


class KPointDomain:
    def __init__(self, k_kc, icell_cv):
        self.k_kc = k_kc
        self.icell_cv = icell_cv

    def __len__(self):
        return len(self.k_kc)

    @cached_property
    def k_kv(self):
        return self.k_kc @ (2 * np.pi * self.icell_cv)
