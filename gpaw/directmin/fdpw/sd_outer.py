"""
Optimization methods for calculating
search directions in space of wafe-functions
Examples are Steepest Descent, Conjugate gradients, L-BFGS
"""


import numpy as np
import copy


class SteepestDescent:
    """
    Steepest descent algorithm
    """

    def __init__(self, wfs, dimensions):
        """
        """
        self.iters = 0
        self.n_kps = wfs.kd.nibzkpts
        self.dimensions = dimensions

    def __str__(self):

        return 'Steepest Descent'

    def update_data(self, psi, g, wfs, prec):
        """
        update search direction

        :param psi:
        :param g:
        :param wfs:
        :param prec:
        :return:
        """

        self.apply_prec(wfs, g, prec, 1.0)

        return self.minus(wfs, g)

    def dot(self, psi_1, psi_2, kpt, wfs):
        """
        dot product between two objects pse_1 and psi_2

        :param psi_1:
        :param psi_2:
        :param kpt:
        :param wfs:
        :return:
        """

        # def S(psit_G):
        #     return psit_G

        def dS(a, P_ni):
            """
            apply PAW
            :param a:
            :param P_ni:
            :return:
            """

            return np.dot(P_ni, wfs.setups[a].dO_ii)

        P1_ai = wfs.pt.dict(shape=1)
        P2_ai = wfs.pt.dict(shape=1)

        wfs.pt.integrate(psi_1, P1_ai, kpt.q)
        wfs.pt.integrate(psi_2, P2_ai, kpt.q)

        dot_prod = wfs.integrate(psi_1, psi_2, False)

        paw_dot_prod = 0.0

        for a in P1_ai.keys():
            paw_dot_prod += \
                np.dot(dS(a, P2_ai[a]), P1_ai[a].T.conj())[0][0]

        sum_dot = dot_prod + paw_dot_prod
        sum_dot = wfs.gd.comm.sum(sum_dot)

        return sum_dot

    def dot_2(self, psi_1, psi_2, kpt, wfs):
        """
        dot product between pseudo-parts

        :param psi_1:
        :param psi_2:
        :param kpt:
        :param wfs:
        :return:
        """

        dot_prod = wfs.integrate(psi_1, psi_2, False)
        dot_prod = wfs.gd.comm.sum(dot_prod)

        return dot_prod

    def dot_all_k_and_b(self, x1, x2, wfs):
        """
        dot product between x1 and x2 over all k-points and bands

        :param x1:
        :param x2:
        :param wfs:
        :return:
        """

        if wfs.dtype is complex:
            dot_pr_x1x2 = 0.0j
        else:
            dot_pr_x1x2 = 0.0

        for kpt in wfs.kpt_u:
            k = self.n_kps * kpt.s + kpt.q
            for i in range(self.dimensions[k]):
                dot_pr_x1x2 += self.dot_2(
                    x1[k][i], x2[k][i], kpt, wfs)

        dot_pr_x1x2 = wfs.kd.comm.sum(dot_pr_x1x2)

        return 2.0 * dot_pr_x1x2.real

    def calc_diff(self, x1, x2, wfs, const_0=1.0, const=1.0):
        """
        calculate difference beetwen x1 and x2

        :param x1:
        :param x2:
        :param wfs:
        :param const_0:
        :param const:
        :return:
        """
        y_k = {}
        for kpt in wfs.kpt_u:
            k = self.n_kps * kpt.s + kpt.q
            y_k[k] = \
                const_0 * x1[k] - \
                const * x2[k]

        return y_k

    def minus(self, wfs, x):

        p = {}
        for kpt in wfs.kpt_u:
            p[self.n_kps * kpt.s + kpt.q] = \
                - x[self.n_kps * kpt.s + kpt.q].copy()

        return p

    def multiply(self, x, const=1.0):

        y = {}
        for k in x.keys():
            y[k] = const * x[k]

        return y

    def zeros(self, x):

        y = {}
        for k in x.keys():
            y[k] = np.zeros_like(x[k])

        return y

    def apply_prec(self, wfs, x, prec, const=1.0):
        """
        apply preconditioning to the gradient

        :param wfs:
        :param x:
        :param prec:
        :param const:
        :return:
        """
        if wfs.mode == 'pw':
            deg = (3.0 - wfs.kd.nspins)
            deg *= 2.0
            for kpt in wfs.kpt_u:
                k = self.n_kps * kpt.s + kpt.q
                for i, y in enumerate(x[k]):
                    psit_G = kpt.psit.array[i]
                    ekin = prec.calculate_kinetic_energy(psit_G, kpt)
                    x[k][i] = - const * prec(y, kpt, ekin) / deg

        else:
            deg = (3.0 - wfs.kd.nspins)
            for kpt in wfs.kpt_u:
                k = self.n_kps * kpt.s + kpt.q
                for i, y in enumerate(x[k]):
                    x[k][i] = - const * prec(y, kpt, None) / deg


class FRcg(SteepestDescent):
    """
    The Fletcher-Reeves conj. grad. method
    See Jorge Nocedal and Stephen J. Wright 'Numerical
    Optimization' Second Edition, 2006 (p. 121)
    """

    def __init__(self, wfs, dimensions):

        """
        """
        super().__init__(wfs, dimensions)

    def __str__(self):
        return 'Fletcher-Reeves conjugate gradient method'

    def update_data(self, psi, g_k1, wfs, prec):

        if prec is not None:
            self.apply_prec(wfs, g_k1, prec, 1.0)

        if self.iters == 0:
            self.p_k = self.minus(wfs, g_k1)
            # save the step
            self.g_k = g_k1
            self.iters += 1
            return self.p_k
        else:
            dot_gg_k1 = self.dot_all_k_and_b(g_k1, g_k1, wfs)
            dot_gg_k = self.dot_all_k_and_b(self.g_k, self.g_k, wfs)
            beta_k = dot_gg_k1 / dot_gg_k
            self.p_k = self.calc_diff(g_k1, self.p_k, wfs,
                                      const_0=-1.0,
                                      const=-beta_k)
            # self.p_k = -g_k1 + beta_k * self.p_k
            # save this step
            self.g_k = g_k1
            self.iters += 1

            if self.iters > 10:
                self.iters = 0

            return self.p_k


class LBFGS(SteepestDescent):
    """
    The limited-memory BFGS.
    See Jorge Nocedal and Stephen J. Wright 'Numerical
    Optimization' Second Edition, 2006 (p. 177)
    """

    def __init__(self, wfs, dimensions, memory=1):

        """
        """
        super().__init__(wfs, dimensions)

        self.s_k = {i: None for i in range(memory)}
        self.y_k = {i: None for i in range(memory)}
        self.rho_k = np.zeros(shape=memory)

        self.kp = {}
        self.p = 0
        self.k = 0  # number of calls

        self.m = memory
        self.stable = True

    def __str__(self):

        return 'LBFGS'

    def update_data(self, x_k1, g_k1, wfs, prec):
        """
        update search direction

        :param x_k1:
        :param g_k1:
        :param wfs:
        :param prec:
        :return:
        """

        if prec is not None:
            self.apply_prec(wfs, g_k1, prec, 1.0)

        if self.k == 0:

            self.kp[self.k] = self.p
            self.x_k = x_k1
            self.g_k = g_k1
            self.s_k[self.kp[self.k]] = self.zeros(g_k1)
            self.y_k[self.kp[self.k]] = self.zeros(g_k1)
            self.k += 1
            self.p += 1
            self.kp[self.k] = self.p
            p = self.minus(wfs, g_k1)

            return p

        else:
            if self.p == self.m:
                self.p = 0
                self.kp[self.k] = self.p

            rho_k = self.rho_k
            kp = self.kp
            k = self.k
            m = self.m

            self.s_k[kp[k]] = self.calc_diff(x_k1, self.x_k, wfs)
            self.y_k[kp[k]] = self.calc_diff(g_k1, self.g_k, wfs)

            dot_ys = self.dot_all_k_and_b(self.y_k[kp[k]],
                                          self.s_k[kp[k]], wfs)
            if abs(dot_ys) > 0.0:
                rho_k[kp[k]] = 1.0 / dot_ys
            else:
                rho_k[kp[k]] = 1.0e16 * np.sign(dot_ys)

            if rho_k[kp[k]] < 0.0:
                # raise Exception('y_k^Ts_k is not positive!')
                self.stable = False
                self.__init__(wfs, self.dimensions, self.m)
                # we could call self.update,
                # but we already applied prec to g
                self.kp[self.k] = self.p
                self.x_k = x_k1
                self.g_k = g_k1
                self.s_k[self.kp[self.k]] = self.zeros(g_k1)
                self.y_k[self.kp[self.k]] = self.zeros(g_k1)
                self.k += 1
                self.p += 1
                self.kp[self.k] = self.p
                p = self.multiply(g_k1, -1.0)

                return p

            # q = np.copy(g_k1)
            q = copy.deepcopy(g_k1)

            alpha = np.zeros(np.minimum(k + 1, m))
            j = np.maximum(-1, k - m)

            for i in range(k, j, -1):

                dot_sq = self.dot_all_k_and_b(self.s_k[kp[i]], q, wfs)

                alpha[kp[i]] = rho_k[kp[i]] * dot_sq

                q = self.calc_diff(q, self.y_k[kp[i]],
                                   wfs, const=alpha[kp[i]])

                # q -= alpha[kp[i]] * y_k[kp[i]]

            try:
                t = np.maximum(1, k - m + 1)

                dot_yy = self.dot_all_k_and_b(self.y_k[kp[t]],
                                              self.y_k[kp[t]], wfs)

                r = self.multiply(q, 1.0 / (rho_k[kp[t]] * dot_yy))

                # r = q / (
                #       rho_k[kp[t]] * np.dot(y_k[kp[t]], y_k[kp[t]]))
            except ZeroDivisionError:
                # r = 1.0e12 * q
                r = self.multiply(q, 1.0e16)

            for i in range(np.maximum(0, k - m + 1), k + 1):

                dot_yr = self.dot_all_k_and_b(self.y_k[kp[i]], r, wfs)

                beta = rho_k[kp[i]] * dot_yr

                r = self.calc_diff(r, self.s_k[kp[i]], wfs,
                                   const=(beta - alpha[kp[i]]))

                # r += s_k[kp[i]] * (alpha[kp[i]] - beta)

            # save this step:
            self.x_k = x_k1
            self.g_k = g_k1

            self.k += 1
            self.p += 1

            self.kp[self.k] = self.p

            return self.multiply(r, const=-1.0)


class LBFGS_P(SteepestDescent):
    """
    The limited-memory BFGS.
    See Jorge Nocedal and Stephen J. Wright 'Numerical
    Optimization' Second Edition, 2006 (p. 177)

    used with preconditioning
    """

    def __init__(self, wfs, dimensions, memory=1):

        """
        """
        super().__init__(wfs, dimensions)

        self.s_k = {i: None for i in range(memory)}
        self.y_k = {i: None for i in range(memory)}
        self.rho_k = np.zeros(shape=memory)

        self.kp = {}
        self.p = 0
        self.k = 0  # number of calls

        self.m = memory
        self.stable = True

    def __str__(self):

        return 'LBFGS'

    def update_data(self, psi, g_k1, wfs, prec):
        """
        update search direction

        :param psi:
        :param g_k1:
        :param wfs:
        :param prec:
        :return:
        """
        if self.k == 0:

            self.kp[self.k] = self.p
            self.x_k = psi
            self.g_k = g_k1
            self.s_k[self.kp[self.k]] = self.zeros(g_k1)
            self.y_k[self.kp[self.k]] = self.zeros(g_k1)
            self.k += 1
            self.p += 1
            self.kp[self.k] = self.p
            p = self.minus(wfs, g_k1)
            if prec is not None:
                self.apply_prec(wfs, p, prec, 1.0)
            return p

        else:
            if self.p == self.m:
                self.p = 0
                self.kp[self.k] = self.p

            s_k = self.s_k
            x_k = self.x_k
            y_k = self.y_k
            g_k = self.g_k

            x_k1 = psi

            rho_k = self.rho_k

            kp = self.kp
            k = self.k
            m = self.m

            s_k[kp[k]] = self.calc_diff(x_k1, x_k, wfs)
            y_k[kp[k]] = self.calc_diff(g_k1, g_k, wfs)

            dot_ys = self.dot_all_k_and_b(y_k[kp[k]],
                                          s_k[kp[k]], wfs)
            if abs(dot_ys) > 0.0:
                rho_k[kp[k]] = 1.0 / dot_ys
            else:
                rho_k[kp[k]] = 1.0e16 * np.sign(dot_ys)

            if rho_k[kp[k]] < 0.0:
                # raise Exception('y_k^Ts_k is not positive!')
                self.stable = False
                self.__init__(wfs, self.dimensions, self.m)
                # we could call self.update,
                # but we already applied prec to g
                self.kp[self.k] = self.p
                self.x_k = x_k1
                self.g_k = g_k1
                self.s_k[self.kp[self.k]] = self.zeros(g_k1)
                self.y_k[self.kp[self.k]] = self.zeros(g_k1)
                self.k += 1
                self.p += 1
                self.kp[self.k] = self.p
                p = self.multiply(g_k1, -1.0)
                if prec is not None:
                    self.apply_prec(wfs, p, prec, 1.0)
                return p

            # q = np.copy(g_k1)
            q = copy.deepcopy(g_k1)

            alpha = np.zeros(np.minimum(k + 1, m))
            j = np.maximum(-1, k - m)

            for i in range(k, j, -1):

                dot_sq = self.dot_all_k_and_b(s_k[kp[i]], q, wfs)

                alpha[kp[i]] = rho_k[kp[i]] * dot_sq

                q = self.calc_diff(q, y_k[kp[i]],
                                   wfs, const=alpha[kp[i]])

                # q -= alpha[kp[i]] * y_k[kp[i]]

            if prec is not None:
                self.apply_prec(wfs, q, prec, 1.0)
                r = q
            else:
                try:
                    t = np.maximum(1, k - m + 1)

                    dot_yy = self.dot_all_k_and_b(y_k[kp[t]],
                                                  y_k[kp[t]], wfs)

                    r = self.multiply(
                        q, 1.0 / (rho_k[kp[t]] * dot_yy))

                    # r = q / (
                    #  rho_k[kp[t]] * np.dot(y_k[kp[t]], y_k[kp[t]]))
                except ZeroDivisionError:
                    # r = 1.0e12 * q
                    r = self.multiply(q, 1.0e16)

            for i in range(np.maximum(0, k - m + 1), k + 1):

                dot_yr = self.dot_all_k_and_b(y_k[kp[i]], r, wfs)

                beta = rho_k[kp[i]] * dot_yr

                r = self.calc_diff(r, s_k[kp[i]], wfs,
                                   const=(beta - alpha[kp[i]]))

                # r += s_k[kp[i]] * (alpha[kp[i]] - beta)

            # save this step:
            self.x_k = x_k1
            self.g_k = g_k1

            self.k += 1
            self.p += 1

            self.kp[self.k] = self.p

            return self.multiply(r, const=-1.0)


class LSR1P(SteepestDescent):
    """
    Limited memory symmetric rank one quasi-newton algorithm.
    arXiv:2006.15922 [physics.chem-ph]
    """

    def __init__(self, wfs, dimensions, memory=10,
                 method='LSR1', phi=None):
        """

        :param wfs:
        :param dimensions:
        :param memory: number of previous steps to store
        :param method:
        :param phi:
        """
        super().__init__(wfs, dimensions)

        self.u_k = {i: None for i in range(memory)}
        self.j_k = {i: None for i in range(memory)}
        self.yj_k = np.zeros(shape=memory)
        self.method = method
        self.phi = phi

        self.phi_k = np.zeros(shape=memory)
        if self.phi is None:
            assert self.method in ['LSR1', 'LP',
                                   'LBofill', 'Linverse_Bofill'], \
                'Value Error'
            if self.method == 'LP':
                self.phi_k.fill(1.0)
        else:
            self.phi_k.fill(self.phi)

        self.kp = {}
        self.p = 0
        self.k = 0

        self.m = memory

    def __str__(self):

        return 'LSR1P'

    def update_data(self, x_k1, g_k1, wfs, precond=None):

        bg_k1 = copy.deepcopy(g_k1)
        if precond is not None:
            self.apply_prec(wfs, bg_k1, precond, 1.0)

        if self.k == 0:
            self.kp[self.k] = self.p
            self.x_k = copy.deepcopy(x_k1)
            self.g_k = copy.deepcopy(g_k1)
            self.u_k[self.kp[self.k]] = self.zeros(g_k1)
            self.j_k[self.kp[self.k]] = self.zeros(g_k1)
            self.k += 1
            self.p += 1
            self.kp[self.k] = self.p
            p = self.minus(wfs, bg_k1)
            self.iters += 1

            return p

        else:
            if self.p == self.m:
                self.p = 0
                self.kp[self.k] = self.p

            x_k = self.x_k
            g_k = self.g_k
            u_k = self.u_k
            j_k = self.j_k
            yj_k = self.yj_k
            phi_k = self.phi_k

            x_k1 = copy.deepcopy(x_k1)

            kp = self.kp
            k = self.k
            m = self.m

            s_k = self.calc_diff(x_k1, x_k, wfs)
            y_k = self.calc_diff(g_k1, g_k, wfs)
            by_k = copy.deepcopy(y_k)

            if precond is not None:
                self.apply_prec(wfs, by_k, precond, 1.0)

            by_k = self.update_bv(
                wfs, by_k, y_k, u_k, j_k, yj_k,
                phi_k, np.maximum(1, k - m), k)

            j_k[kp[k]] = self.calc_diff(s_k, by_k, wfs)
            yj_k[kp[k]] = self.dot_all_k_and_b(y_k, j_k[kp[k]], wfs)

            dot_yy = self.dot_all_k_and_b(y_k, y_k, wfs)
            if abs(dot_yy) > 1.0e-15:
                u_k[kp[k]] = self.multiply(y_k, 1.0 / dot_yy)
            else:
                u_k[kp[k]] = self.multiply(y_k, 1.0e15)

            if self.method == 'LBofill' and self.phi is None:
                jj_k = self.dot_all_k_and_b(
                    j_k[kp[k]], j_k[kp[k]], wfs)
                phi_k[kp[k]] = 1 - yj_k[kp[k]]**2 / (dot_yy * jj_k)
            elif self.method == 'Linverse_Bofill' and \
                    self.phi is None:
                jj_k = self.dot_all_k_and_b(
                    j_k[kp[k]], j_k[kp[k]], wfs)
                phi_k[kp[k]] = yj_k[kp[k]] ** 2 / (dot_yy * jj_k)

            bg_k1 = self.update_bv(
                wfs, bg_k1, g_k1, u_k, j_k, yj_k, phi_k,
                np.maximum(1, k - m + 1), k + 1)

            # save this step:
            self.x_k = copy.deepcopy(x_k1)
            self.g_k = copy.deepcopy(g_k1)
            self.k += 1
            self.p += 1
            self.kp[self.k] = self.p
            self.iters += 1

        return self.multiply(bg_k1, const=-1.0)

    def update_bv(self, wfs, bv, v, u_k, j_k, yj_k, phi_k, i_0, i_m):
        kp = self.kp

        for i in range(i_0, i_m):
            dot_uv = self.dot_all_k_and_b(u_k[kp[i]],
                                          v, wfs)
            dot_jv = self.dot_all_k_and_b(j_k[kp[i]],
                                          v, wfs)

            alpha = dot_jv - yj_k[kp[i]] * dot_uv
            beta_p = self.calc_diff(j_k[kp[i]], u_k[kp[i]],
                                    wfs, const_0=dot_uv,
                                    const=-alpha)

            beta_ms = self.multiply(j_k[kp[i]], dot_jv / yj_k[kp[i]])

            beta = self.calc_diff(beta_ms, beta_p, wfs,
                                  const_0=1.0 - phi_k[kp[i]],
                                  const=-phi_k[kp[i]])

            bv = self.calc_diff(bv, beta, wfs, const=-1.0)

        return bv
