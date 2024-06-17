import numpy as np
from functools import cached_property


class ResponseKPointGrid:
    def __init__(self, kd, icell_cv, bzk_kc=None):
        # Very hacky that bzk_kc may or may not be the same as those on kd.
        # Sometimes the code likes to "just have an array", maybe to avoid
        # depending on whatever arbitrary processing is hardcoded on kd.

        self.icell_cv = icell_cv
        if bzk_kc is None:
            bzk_kc = kd.bzk_kc
        self.bzk_kc = bzk_kc
        self.bzk_kv = bzk_kc @ (2 * np.pi * icell_cv)

        # XXX May or may not be consistent with bzk_kc!
        self.kd = kd

    @cached_property
    def kptfinder(self):
        from gpaw.response.symmetry import KPointFinder
        return KPointFinder(self.bzk_kc)


class KPointDomain:
    def __init__(self, k_kc, icell_cv):
        self.k_kc = k_kc
        self.icell_cv = icell_cv

    def __len__(self):
        return len(self.k_kc)

    @cached_property
    def k_kv(self):
        return self.k_kc @ (2 * np.pi * self.icell_cv)
