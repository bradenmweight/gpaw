#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from math import pi, log, sqrt

import numpy as np
from scipy.special import gamma
from scipy.interpolate import interp1d

from ase.units import Hartree
from ase.utils import prnt

from gpaw.atom.configurations import configurations
from gpaw.atom.aeatom import AllElectronAtom, Channel, parse_ld_str, colors
from gpaw.setup import BaseSetup
from gpaw.spline import Spline
from gpaw.basis_data import Basis
#from gpaw.hgh import null_xc_correction
from gpaw.setup_data import SetupData
from gpaw.version import version


class PAWWaves:
    def __init__(self, rgd, l, rcut):
        self.rgd = rgd
        self.l = l
        self.rcut = rcut

        self.n_n = []
        self.e_n = []
        self.f_n = []
        self.phi_ng = []
        self.phit_ng = None
        self.pt_ng = None

    def __len__(self):
        return len(self.n_n)
    
    def add(self, phi_g, n, e, f):
        self.phi_ng.append(phi_g)
        self.n_n.append(n)
        self.e_n.append(e)
        self.f_n.append(f)

    def pseudize(self):
        rgd = self.rgd

        phi_ng = self.phi_ng = np.array(self.phi_ng)
        N = len(phi_ng)
        phit_ng = self.phit_ng = rgd.empty(N)
        gcut = rgd.ceil(self.rcut)

        self.nt_g = 0
        self.c_np = []
        for n in range(N):
            phit_ng[n], c_p = rgd.pseudize(phi_ng[n], gcut, self.l, points=6)
            self.c_np.append(c_p)
            self.nt_g += self.f_n[n] / 4 / pi * phit_ng[n]**2
            
        self.dS_nn = np.empty((N, N))
        for n1 in range(N):
            for n2 in range(N):
                self.dS_nn[n1, n2] = rgd.integrate(
                    phi_ng[n1] * phi_ng[n2] -
                    phit_ng[n1] * phit_ng[n2]) / (4 * pi)
        self.Q = np.dot(self.f_n, self.dS_nn.diagonal())

    def construct_projectors(self, vtr_g):
        rgd = self.rgd
        phit_ng = self.phit_ng
        N = len(phit_ng)
        gcut = rgd.ceil(self.rcut)
        r_g = rgd.r_g
        l = self.l
        P = len(self.c_np[0]) - 1
        p = np.arange(2 * P, 0, -2) + l
        A_nn = np.empty((N, N))
        q_ng = rgd.zeros(N)
        for n in range(N):
            q_g = ((vtr_g - self.e_n[n] * r_g) * self.phit_ng[n] +
                   np.polyval(-0.5 * self.c_np[n][:P] *
                               (p * (p + 1) - l * (l + 1)), r_g**2) *
                   r_g**(1 + l))
            q_g[gcut:] = 0
            q_ng[n] = q_g
            A_nn[n] = rgd.integrate(q_g * phit_ng, -1) / (4 * pi)

        self.pt_ng = np.dot(np.linalg.inv(A_nn), q_ng)
        self.pt_ng[:, 1:] /= r_g[1:]
        self.pt_ng[:, 0] = self.pt_ng[:, 1]
        self.dH_nn = self.e_n * self.dS_nn - A_nn.T

    def calculate_kinetic_energy_correction(self, vr_g, vtr_g):
        self.dekin_nn = (self.rgd.integrate(self.phi_ng[:, np.newaxis] *
                                           self.phi_ng *
                                           vr_g, -1) / (4 * pi) -
                         self.rgd.integrate(self.phit_ng[:, np.newaxis] *
                                           self.phit_ng *
                                           vtr_g, -1) / (4 * pi) -
                        self.dH_nn)


class PAWSetupGenerator:
    def __init__(self, aea, states, radii, fd=sys.stdout):
        """fd: stream
            Text output."""
        
        self.aea = aea
        self.rgd = aea.rgd
        
        if fd is None:
            fd = devnull
        self.fd = fd

        self.log('Generating PAW setup.')

        self.rcmax = max(radii)
        self.gcmax = self.rgd.ceil(self.rcmax)

        self.waves_l = []
        for l, nn in enumerate(states):
            waves = PAWWaves(self.rgd, l, radii[l])
            for n in nn:
                if isinstance(n, int):
                    # Bound state:
                    ch = aea.channels[l]
                    e = ch.e_n[n - l - 1]
                    f = ch.f_n[n - l - 1]
                    phi_g = ch.phi_ng[n - l - 1]
                else:
                    e = n
                    n = -1
                    f = 0.0
                    phi_g = self.rgd.zeros()
                    gc = self.gcmax + 10
                    ch = Channel(l)
                    ch.integrate_outwards(phi_g, self.rgd, aea.vr_sg[0], gc, e)
                    phi_g[1:gc] /= self.rgd.r_g[1:gc]
                    if l == 0:
                        phi_g[0] = phi_g[1]
                waves.add(phi_g, n, e, f)
            self.waves_l.append(waves)

        self.alpha = log(1.0e4) / self.rcmax**2  # exp(-alpha*rcmax^2)=1.0e-4
        self.alpha = round(self.alpha, 2)
        self.ghat_g = (np.exp(-self.alpha * self.rgd.r_g**2) *
                       (self.alpha / pi)**1.5)
        
        self.vtr_g = None

    def log(self, *args, **kwargs):
        prnt(file=self.fd, *args, **kwargs)

    def calculate_core_density(self):
        self.nc_g = self.rgd.zeros()
        self.ncore = 0
        self.nvalence = 0
        self.ekincore = 0.0
        for l, ch in enumerate(self.aea.channels):
            for n, f in enumerate(ch.f_n):
                if (l >= len(self.waves_l) or
                    (l < len(self.waves_l) and
                     n + l + 1 not in self.waves_l[l].n_n)):
                    self.nc_g += f * ch.calculate_density(n)
                    self.ncore += f
                    self.ekincore += f * ch.e_n[n]
                else:
                    self.nvalence += f
        
        self.ekincore -= self.rgd.integrate(self.nc_g * self.aea.vr_sg[0], -1)
        self.nct_g = self.rgd.pseudize(self.nc_g, self.gcmax)[0]
        self.npseudocore = self.rgd.integrate(self.nct_g)

        self.log('Core electrons:', self.ncore)
        self.log('Pseudo core electrons: %.3f' % self.rgd.integrate(self.nct_g))
        
        self.nt_g = self.nct_g.copy()
        self.Q = -self.aea.Z + self.ncore - self.npseudocore

    def pseudize(self):
        for waves in self.waves_l:
            waves.pseudize()
            self.nt_g += waves.nt_g
            self.Q += waves.Q

        self.rhot_g = self.nt_g + self.Q * self.ghat_g
        self.vHtr_g = self.rgd.poisson(self.rhot_g)

        self.vxct_g = self.rgd.zeros()
        exct_g = self.rgd.zeros()
        self.exct = self.aea.xc.calculate_spherical(
            self.rgd, self.nt_g.reshape((1, -1)), self.vxct_g.reshape((1, -1)))

        self.log('\nProjectors:')
        self.log(' state  occ         energy             norm        rcut')
        self.log(' nl            [Hartree]  [eV]      [electrons]   [Bohr]')
        self.log('----------------------------------------------------------')
        for l, waves in enumerate(self.waves_l):
            for n, e, f, ds in zip(waves.n_n, waves.e_n, waves.f_n,
                                  waves.dS_nn.diagonal()):
                if n == -1:
                    self.log('  %s             %10.6f %10.5f %10.2f' %
                             ('spdf'[l], e, e * Hartree, waves.rcut))
                else:
                    self.log(
                        ' %d%s     %2d       %10.6f %10.5f    %5.3f %9.2f' %
                             (n, 'spdf'[l], f, e, e * Hartree, 1 - ds,
                              waves.rcut))
                    
    def find_polynomial_potential(self, r0, P, e0=None):
        g0 = self.rgd.ceil(r0)
        vtr_g = self.vHtr_g + self.vxct_g * self.rgd.r_g
        if e0 is None:
            self.vtr_g = self.rgd.pseudize(vtr_g, g0, 1, P)[0]
            self.v0r_g = self.vtr_g - vtr_g
        else:
            asdgfkljh

    def find_local_potential(self, l0, r0, P, e0):
        self.log('Local potential matching %s-state at %.3f eV' %
                 ('spdfg'[l0], e0 * Hartree))
        
        g0 = self.rgd.ceil(r0)
        gc = g0 + 20

        ch = Channel(l0)
        phi_g = self.rgd.zeros()
        ch.integrate_outwards(phi_g, self.rgd, self.aea.vr_sg[0], gc, e0)
        phi_g[1:gc] /= self.rgd.r_g[1:gc]
        if l0 == 0:
            phi_g[0] = phi_g[1]

        phit_g, c_p = self.rgd.pseudize_normalized(phi_g, g0, l=l0, points=P)
        r_g = self.rgd.r_g[1:g0]
        p = np.arange(2 * P, 0, -2) + l0
        t_g = np.polyval(-0.5 * c_p[:P] * (p * (p + 1) - l0 * (l0 + 1)),
                          r_g**2) 
        self.vtr_g = self.aea.vr_sg[0].copy()
        self.vtr_g[0] = 0.0
        self.vtr_g[1:g0] = e0 * r_g - t_g * r_g**(l0 + 1) / phit_g[1:g0]
        self.v0r_g = self.vtr_g - self.vHtr_g - self.vxct_g * self.rgd.r_g

    def construct_projectors(self):
        for waves in self.waves_l:
            waves.construct_projectors(self.vtr_g)
            waves.calculate_kinetic_energy_correction(self.aea.vr_sg[0],
                                                      self.vtr_g)

    def plot(self):
        import matplotlib.pyplot as plt
        r_g = self.rgd.r_g

        plt.figure()
        plt.plot(r_g, self.vxct_g, label='xc')
        plt.plot(r_g[1:], self.v0r_g[1:] / r_g[1:], label='0')
        plt.plot(r_g[1:], self.vHtr_g[1:] / r_g[1:], label='H')
        plt.plot(r_g[1:], self.vtr_g[1:] / r_g[1:], label='ps')
        plt.plot(r_g[1:], self.aea.vr_sg[0, 1:] / r_g[1:], label='ae')
        plt.axis(xmax=2 * self.rcmax,
                 ymin=self.vtr_g[1] / r_g[1],
                 ymax=max(0, (self.v0r_g[1:] / r_g[1:]).max()))
        plt.legend()
        
        plt.figure()
        for waves in self.waves_l:
            for phi_g, phit_g in zip(waves.phi_ng, waves.phit_ng):
                plt.plot(r_g, phi_g * r_g)
                plt.plot(r_g, phit_g * r_g, '--')
        plt.axis(xmax=3 * self.rcmax)

        plt.figure()
        for waves in self.waves_l:
            for pt_g in waves.pt_ng:
                plt.plot(r_g, pt_g * r_g)
        plt.axis(xmax=3 * self.rcmax)

    def logarithmic_derivative(self, l, energies, rcut):
        rgd = self.rgd
        ch = Channel(l)
        gcut = rgd.round(rcut)

        if l < len(self.waves_l):
            # Nonlocal PAW stuff:
            waves = self.waves_l[l]
            pt_ng = waves.pt_ng
            dH_nn = waves.dH_nn
            dS_nn = waves.dS_nn
            N = len(pt_ng)
        else:
            N = 0

        u_g = rgd.zeros()
        u_ng = rgd.zeros(N)
        dudr_n = np.empty(N)

        logderivs = []
        for e in energies:
            dudr = ch.integrate_outwards(u_g, rgd, self.vtr_g, gcut, e)
            u = u_g[gcut]
            
            if N:
                for n in range(N):
                    dudr_n[n] = ch.integrate_outwards(u_ng[n], rgd,
                                                      self.vtr_g, gcut, e,
                                                      pt_ng[n])
            
                A_nn = dH_nn - e * dS_nn
                B_nn = rgd.integrate(pt_ng * u_ng[:, np.newaxis], -1)
                c_n  = rgd.integrate(pt_ng * u_g, -1)
                d_n = np.linalg.solve(np.dot(A_nn, B_nn) - 4 * pi * np.eye(N),
                                      np.dot(A_nn, c_n))
                u -= np.dot(u_ng[:, gcut], d_n)
                dudr -= np.dot(dudr_n, d_n)
            
            logderivs.append(dudr / u)

        return logderivs
    
    def make_paw_setup(self):
        aea = self.aea
        
        setup = SetupData(aea.symbol, aea.xc.name, 'paw', readxml=False)

        nj = sum(len(waves) for waves in self.waves_l)
        setup.e_kin_jj = np.zeros((nj, nj))
        setup.id_j = []
        j1 = 0
        for l, waves in enumerate(self.waves_l):
            ne = 0
            for n, f, e, phi_g, phit_g, pt_g in zip(waves.n_n, waves.f_n,
                                                    waves.e_n, waves.phi_ng,
                                                    waves.phit_ng,
                                                    waves.pt_ng):
                setup.append(n, l, f, e, waves.rcut, phi_g, phit_g, pt_g)
                if n == -1:
                    ne += 1
                    id = '%s-%s%d' % (aea.symbol, 'spdf'[l], ne)
                else:
                    id = '%s-%d%s' % (aea.symbol, n, 'spdf'[l])
                setup.id_j.append(id)
            j2 = j1 + len(waves)
            setup.e_kin_jj[j1:j2, j1:j2] = waves.dekin_nn
            j1 = j2

        setup.nc_g = self.nc_g * sqrt(4 * pi)
        setup.nct_g = self.nct_g * sqrt(4 * pi)
        setup.e_kinetic_core = self.ekincore
        setup.vbar_g = self.v0r_g * sqrt(4 * pi)
        setup.vbar_g[1:] /= self.rgd.r_g[1:]
        setup.vbar_g[0] = setup.vbar_g[1]
        setup.Z = aea.Z
        setup.Nc = self.ncore
        setup.Nv = self.nvalence
        setup.e_kinetic = aea.ekin
        setup.e_xc = aea.exc
        setup.e_electrostatic = aea.eH + aea.eZ
        setup.e_total = aea.exc + aea.ekin + aea.eH + aea.eZ
        setup.rgd = self.rgd
        setup.rcgauss = 1 / sqrt(self.alpha)

        setup.tauc_g = self.rgd.zeros()
        setup.tauct_g = self.rgd.zeros()
        print 'no tau!!!!!!!!!!!'
        
        reltype = 'non-relativistic'
        attrs = [('type', reltype), ('name', 'gpaw-%s' % version)]
        setup.generatorattrs = attrs

        return setup

     
def build_parser(): 
    from optparse import OptionParser

    parser = OptionParser(usage='%prog [options] element',
                          version='%prog 0.1')
    parser.add_option('-f', '--xc-functional', type='string', default='LDA',
                      help='Exchange-Correlation functional ' +
                      '(default value LDA)',
                      metavar='<XC>')
    parser.add_option('-P', '--projectors',
                      help='Projector functions - use comma-separated - ' +
                      'nl values, where n can be pricipal quantum number ' +
                      '(integer) or energy (floating point number). ' +
                      'Example: 2s,0.5s,2p,0.5p,0.0d.')
    parser.add_option('-r', '--radius',
                      help='1.2 or 1.2,1.1,1.1')
    parser.add_option('-0', '--zero-potential',
                      metavar='type,nderivs,radius,e0',
                      help='Parameters for zero potential.')
    parser.add_option('-p', '--plot', action='store_true')
    parser.add_option('-l', '--logarithmic-derivatives',
                      metavar='spdfg,e1:e2:de,radius',
                      help='Plot logarithmic derivatives. ' +
                      'Example: -l spdf,-1:1:0.05,1.3. ' +
                      'Energy range and/or radius can be left out.')
    parser.add_option('-w', '--write', action='store_true')
    return parser


def main(AEA=AllElectronAtom):
    parser = build_parser()
    opt, args = parser.parse_args()

    if len(args) != 1:
        parser.error('Incorrect number of arguments')
    symbol = args[0]

    kwargs = {'xc': opt.xc_functional}
        
    aea = AEA(symbol, **kwargs)

    lmax = -1
    states = {}
    if opt.projectors:
        for s in opt.projectors.split(','):
            l = 'spdf'.find(s[-1])
            n = eval(s[:-1])
            if l in states:
                states[l].append(n)
            else:
                states[l] = [n]
            if l > lmax:
                lmax = l

    # Add empty bound states:
    for l, nn in states.items():
        for n in nn:
            if (isinstance(n, int) and
                (l not in aea.f_lsn or n - l > len(aea.f_lsn[l][0]))):
                aea.add(n, l, 0)

    aea.initialize()
    aea.run()
    aea.refine()
                
    if opt.radius:
        radii = [float(r) for r in opt.radius.split(',')]
    else:
        rave = aea.rgd.integrate(aea.n_sg[0], 1) / aea.Z
        radii = [0.5 * rave]
    if lmax >= 0:
        radii += [radii[-1]] * (lmax + 1 - len(radii))

    projectors = [states[l] for l in range(lmax + 1)]

    gen = PAWSetupGenerator(aea, projectors, radii)

    gen.calculate_core_density()
    gen.pseudize()

    if opt.zero_potential:
        x = opt.zero_potential.split(',')
        type = x[0]
        nderiv = int(x[1])
        r0 = float(x[2])
        if len(x) == 4:
            e0 = float(x[3])
        elif type == 'poly':
            e0 = None
        else:
            e0 = 0.0

        if type == 'poly':
            gen.find_polynomial_potential(r0, nderiv, e0)
        else:
            l0 = 'spdfg'.find(type)
            gen.find_local_potential(l0, r0, nderiv, e0)
    else:
        gen.find_polynomial_potential(max(radii), 4)
        
    gen.construct_projectors()

    if opt.write:
        gen.make_paw_setup().write_xml()
        
    if opt.logarithmic_derivatives or opt.plot:
        import matplotlib.pyplot as plt
        if opt.logarithmic_derivatives:
            r = 1.1 * max(radii)
            lvalues, energies, r = parse_ld_str(opt.logarithmic_derivatives, r)
            for l in lvalues:
                ld = aea.logarithmic_derivative(l, energies, r)
                plt.plot(energies, ld, colors[l], label='spdfg'[l])
                ld = gen.logarithmic_derivative(l, energies, r)
                plt.plot(energies, ld, '--' + colors[l])
            plt.legend(loc='best')

        if opt.plot:
            gen.plot()

        plt.show()


if __name__ == '__main__':
    main()

class PAWSetup:
    def __init__(self, alpha, r_g, phit_g, v0_g):
        self.natoms = 0
        self.E = 0.0
        self.Z = 1
        self.Nc = 0
        self.Nv = 1
        self.niAO = 1
        self.pt_j = []
        self.ni = 0
        self.l_j = []
        self.nct = None
        self.Nct = 0.0

        rc = 1.0
        r2_g = np.linspace(0, rc, 100)**2
        x_g = np.exp(-alpha * r2_g)
        x_g[-1] = 0 

        self.ghat_l = [Spline(0, rc,
                              (4 * pi)**0.5 * (alpha / pi)**1.5 * x_g)]

        self.vbar = Spline(0, rc, (4 * pi)**0.5 * v0_g[0] * x_g)

        r = np.linspace(0, 4.0, 100)
        phit = splev(r, splrep(r_g, phit_g))
        poly = np.polyfit(r[[-30,-29,-2,-1]], [0, 0, phit[-2], phit[-1]], 3)
        phit[-30:] -= np.polyval(poly, r[-30:])
        self.phit_j = [Spline(0, 4.0, phit)]
                              
        self.Delta_pL = np.zeros((0, 1))
        self.Delta0 = -1 / (4 * pi)**0.5
        self.lmax = 0
        self.K_p = self.M_p = self.MB_p = np.zeros(0)
        self.M_pp = np.zeros((0, 0))
        self.Kc = 0.0
        self.MB = 0.0
        self.M = 0.0
        self.xc_correction = null_xc_correction
        self.HubU = None
        self.dO_ii = np.zeros((0, 0))
        self.type = 'local'
        self.fingerprint = None
        
    def get_basis_description(self):
        return '1s basis cut off at 4 Bohr'

    def print_info(self, text):
        text('Local pseudo potential')
        
    def calculate_initial_occupation_numbers(self, magmom, hund, charge,
                                             nspins):
        return np.array([(1.0,)])

    def initialize_density_matrix(self, f_si):
        return np.zeros((1, 0))

    def calculate_rotations(self, R_slmm):
        self.R_sii = np.zeros((1, 0, 0))

if 0:
    aea = AllElectronAtom('H')
    aea.add(2, 1, 0)
    aea.initialize()
    aea.run()
    g = PAWSetupGenerator(aea, [(0.8, [1])], 1)
    #g = PAWSetupGenerator(aea, [], 0, 0.8)
    g.generate()
    g.plot()
    setup = g.make_paw_setup()
    from ase.data.molecules import molecule
    from gpaw import GPAW

    a = molecule('H', pbc=1, magmoms=[0])
    a.center(vacuum=2)
    a.set_calculator(
        GPAW(setups={0: setup}))
    a.get_potential_energy()
