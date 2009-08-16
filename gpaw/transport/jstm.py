import time
import numpy as np
import numpy.linalg as la
from ase.units import Bohr, Hartree
from gpaw.operators import Laplace
from gpaw.grid_descriptor import GridDescriptor
from gpaw.lfc import NewLocalizedFunctionsCollection as LFC
from gpaw.utilities.tools import tri2full
import gpaw.mpi as mpi
from gpaw.lcao.tools import remove_pbc, get_lead_lcao_hamiltonian
from math import cos, sin, pi
import cPickle as pickle
from gpaw.lcao.projected_wannier import dots
from gpaw import GPAW
from gpaw.mpi import world as w

class LocalizedFunctions:
    def __init__(self, gd, f_iG, corner_c, index=None, vt_G=None):
        self.gd = gd
        self.size_c = np.array(f_iG.shape[1:4])
        self.f_iG = f_iG
        self.corner_c = corner_c
        self.index = index
        self.vt_G = vt_G
        self.restricted = False
        self.phase_vector_c = np.array([0,0])

    def __len__(self):
        return len(self.f_iG)

    def set_phase_vector(self, v_c):
        self.phase_vector_c = np.asarray(v_c)

    def apply_t(self):
        """Apply kinetic energy operator and return new object."""
        p = 2  # padding
        newsize_c = self.size_c + 2 * p
        gd = GridDescriptor(N_c=newsize_c + 1,
                            cell_cv=self.gd.h_c * (newsize_c + 1),
                            pbc_c=False,
                            comm=mpi.serial_comm)
        T = Laplace(gd, scale=-1/2., n=p)
        f_ig = np.zeros((len(self.f_iG),) + tuple(newsize_c))
        f_ig[:, p:-p, p:-p, p:-p] = self.f_iG
        Tf_iG = np.empty_like(f_ig)
        T.apply(f_ig, Tf_iG)
        return LocalizedFunctions(self.gd, Tf_iG, self.corner_c - p,
                                  self.index)
        
    def overlap(self, other):
        start_c = np.maximum(self.corner_c, other.corner_c)
        stop_c = np.minimum(self.corner_c + self.size_c,
                            other.corner_c + other.size_c)
        if (start_c < stop_c).all():
            astart_c = start_c - self.corner_c
            astop_c = stop_c - self.corner_c
            a_iG = self.f_iG[:,
                astart_c[0]:astop_c[0],
                astart_c[1]:astop_c[1],
                astart_c[2]:astop_c[2]].reshape((len(self.f_iG), -1))
            bstart_c = start_c - other.corner_c
            bstop_c = stop_c - other.corner_c
            b_iG = other.f_iG[:,
                bstart_c[0]:bstop_c[0],
                bstart_c[1]:bstop_c[1],
                bstart_c[2]:bstop_c[2]].reshape((len(other.f_iG), -1))
            a_iG1 = a_iG.copy()
            if self.vt_G is not None:
                a_iG1 *= self.vt_G[start_c[0]:stop_c[0],
                                  start_c[1]:stop_c[1],
                                  start_c[2]:stop_c[2]].reshape((-1,))
            return self.gd.dv * np.inner(a_iG1, b_iG)
        else:
            return None
        
    def restrict(self):
        """Restricts the box of the objet to the current grid"""
        start_c = np.maximum(self.corner_c, np.zeros(3))
        stop_c = np.minimum(self.corner_c + self.size_c, self.gd.N_c)
        if (start_c < stop_c).all():
            astart_c = start_c - self.corner_c
            astop_c = stop_c -self.corner_c
            a_iG = self.f_iG[:,
                astart_c[0]:astop_c[0],
                astart_c[1]:astop_c[1],
                astart_c[2]:astop_c[2]]
            new_corner_c = self.corner_c
            for i in np.where(self.corner_c<0):
                new_corner_c[i] = 0
            
            if (self.f_iG.shape != a_iG.shape):
                self.restricted = True
 
            self.corner_c = new_corner_c           
            self.f_iG = a_iG
            self.size_c=np.asarray(a_iG.shape[1:])

    def __or__(self, other):
        if isinstance(other, LocalizedFunctions):
            return self.overlap(other)

        # other is a potential:
        vt_G = other
        return LocalizedFunctions(self.gd, self.f_iG, self.corner_c,
                                  self.index, vt_G)

class WannierFunction(LocalizedFunctions):
    def __init__(self, gd, wanf_G, corner_c, index=None):
        LocalizedFunctions.__init__(self, gd, wanf_G[np.newaxis, :, :, :],
                                    corner_c, index)

class AtomCenteredFunctions(LocalizedFunctions):
    def __init__(self, gd, spline_j, spos_c, index=None):
        rcut = max([spline.get_cutoff() for spline in spline_j])
        self.spos_c = spos_c
        # Determination of vector distance, c, from position to origin of the box
        cell = gd.cell_cv.copy()
        diagonal = (cell[0,:]+cell[1,:])
        diagonal = diagonal/np.linalg.norm(diagonal)
        a = np.zeros_like(diagonal)
        a[0]=cell[0,:][1]
        a[1]=-cell[0,:][0]
        a=-a/np.linalg.norm(a)
        c = rcut/np.dot(diagonal,a)
        # Determine corner
        A = cell.T/gd.cell_c # Basis change matrix
        pos =np.dot(np.linalg.inv(cell.T), diagonal * c)
        pos[2] = rcut/gd.cell_c[2]
        corner_c = np.ceil(spos_c * gd.N_c - pos*gd.cell_c / gd.h_c).astype(int)
        self.center = np.ceil(pos*gd.cell_c / gd.h_c).astype(int)
        size_c = np.ceil(spos_c * gd.N_c + \
                 pos*gd.cell_c / gd.h_c).astype(int) - corner_c

        smallgd = GridDescriptor(N_c=size_c + 1,
                                 cell_cv=(np.dot(A,np.diag(gd.h_c * (size_c +1))).T),
                                 pbc_c=False,
                                 comm=mpi.serial_comm)

        self.test = (np.dot(A,np.diag(gd.h_c * (size_c +1))).T)
        self.smallgd=smallgd
        lfc = LFC(smallgd, [spline_j])
        lfc.set_positions((spos_c[np.newaxis, :] * gd.N_c - corner_c + 1) /
                          smallgd.N_c)
        ni = lfc.Mmax
        f_iG = smallgd.zeros(ni)
        lfc.add(f_iG, {0: np.eye(ni)})
        LocalizedFunctions.__init__(self, gd, f_iG, corner_c, 
                                     index=index)


class STM:
    def __init__(self, tip, surface, lead1 = None, lead2 = None, **kwargs):
        self.tip = tip
        self.srf = surface
        self.lead1 = lead1
        self.lead2 = lead2
        self.stm_calc = None
        self.ediff = None        
        self.scans ={}
        tgd = tip.gd
        sgd = surface.gd

        self.input_parameters = {'tip_atom_index': 0,
                                 'dmin': 6.0,
                                 'hs1': None,
                                 'hs10': None,
                                 'hs2': None,
                                 'hs20': None,
                                 'align_bf': 0,
                                 'cvl1': 0,
                                 'cvl2': 0,
                                 'bias': 1.0,
                                 'de': 0.01,
                                 'energies': None,
                                 'w': 0.0,
                                 'eta1': 1e-3,
                                 'eta2': 1e-3,
                                 'logfile': '-',
                                 'verbose': False}

        self.initialized = False
        self.transport_uptodate = False
        self.set(**kwargs)

    def set(self, **kwargs):
        for key in kwargs:
            if key in ['hs1', 'hs10', 'hs2', 'hs20',
                       'cvl1', 'cvl2', 'bias',
                       'de', 'energies', 'w',
                       'align_bf','eta1','eta2']:
                self.transport_uptodate = False
                break
            elif key in ['tip_atom_index', 'dmin']:
                self.initialized = False
            elif key not in self.input_parameters:
                raise KeyError, '\'%s\' not a valid keyword' % key

        self.input_parameters.update(kwargs)
        log = self.input_parameters['logfile']
        if log is None:
            class Trash:
                def write(self,s):
                    pass
                def flush(self):
                    pass
            self.log = Trash()
        elif log == '-':
            from sys import stdout
            self.log = stdout
        elif 'logfile' in kwargs:
            self.log = open(log, 'w')

    def initialize(self):
        if self.initialized and self.transport_uptodate:
            return
        elif not self.transport_uptodate and self.initialized:
            self.initialize_transport()
            return

        T = time.localtime()
        print >> self.log, '#%d:%02d:%02d'\
                     % (T[3], T[4], T[5])+' Initializing'

        p = self.input_parameters        
        self.dmin = p['dmin'] / Bohr
        tip_atom_index = p['tip_atom_index']   
        
        # preselect tip and surface functions
        tip_pos_av = self.tip.atoms.get_positions() / Bohr
        srf_pos_av = self.srf.atoms.get_positions() / Bohr
        tip_zmin = tip_pos_av[tip_atom_index, 2]
        srf_zmax = srf_pos_av[:, 2].max()
        tip_zmin_a = np.empty(len(tip_pos_av))
        
        for a, setup in enumerate(self.tip.wfs.setups):
            rcutmax = max([phit.get_cutoff() for phit in setup.phit_j])
            tip_zmin_a[a] = tip_pos_av[a, 2] - rcutmax - tip_zmin
        srf_zmax_a = np.empty(len(srf_pos_av))
        for a, setup in enumerate(self.srf.wfs.setups):
            rcutmax = max([phit.get_cutoff() for phit in setup.phit_j])
            srf_zmax_a[a] = srf_pos_av[a, 2] + rcutmax - srf_zmax
        
        tip_indices = np.where(tip_zmin_a < srf_zmax_a.max() - self.dmin)[0]  
        srf_indices = np.where(srf_zmax_a > tip_zmin_a.min() + self.dmin)[0]  
 
        print >> self.log, '# dmin = ', str(self.dmin*Bohr)
        print >> self.log, '#Tip atoms: ', str(tip_indices.min()) + ' to ' + str(tip_indices.max())
        print >> self.log, '#Surface atoms: ', str(srf_indices.min()) + ' to ' + str(srf_indices.max())

        # tip initialization
        self.tip_cell = TipCell(self.tip, self.srf)
        self.tip_cell.initialize(tip_indices, tip_atom_index)
        self.ni = self.tip_cell.ni       

        # surface initialization
        self.srf_cell = SrfCell(self.srf)
        self.srf_cell.initialize(self.tip_cell, srf_indices)
        self.nj = self.srf_cell.nj
        
        if not self.transport_uptodate:
            self.initialize_transport()            

        self.set_tip_position([0, 0])
        self.initialized = True
        self.log.flush()

    def initialize_transport(self, restart= False, ediff_only = False):
        p = self.input_parameters        
        h1, s1 = p['hs1']
        h10, s10 = p['hs10']
        h2, s2 = p['hs2']
        h20, s20 = p['hs20']
        cvl1 = p['cvl1']
        cvl2 = p['cvl2']
        align_bf = p['align_bf']
        de = p['de']
        bias = p['bias']        
        w = p['w']
        eta1 = p['eta1']
        eta2 = p['eta2']
        
        if cvl1 == 0:
            cvl1 = 1
            
        h1 = h1[:-cvl1, :-cvl1]
        s1 = s1[:-cvl1, :-cvl1]
        h2 = h2[cvl2:, cvl2:]
        s2 = s2[cvl2:, cvl2:]

        # Align bfs with the surface lead as a reference
        diff = (h2[align_bf,align_bf] - h20[align_bf, align_bf]) \
               / s2[align_bf, align_bf]
        diff = diff.real
        h2 -= diff * s2      
        h1 -= diff * s1        
        
        self.ediff = diff + (self.tip.get_fermi_level()\
                + self.srf.get_fermi_level())/2.0

        diff1 = (h10[-1,-1] - h1[-1,-1])/s1[-1,-1]
        h10 -= diff1 * s10
        if not ediff_only and not self.transport_uptodate:
            from ase.transport.stm import STM as STMCalc
            T = time.localtime()
            print >> self.log, '#%d:%02d:%02d'\
                     % (T[3], T[4], T[5])+' Updating transport calculator'

            if p['energies'] == None:
                energies = np.sign(bias) * \
                np.arange(-abs(bias)*w, -abs(bias)*(w-1)+de, de)
                energies.sort()
            else:
                energies = p['energies']

            stm_calc = STMCalc(h2,  s2, 
                               h1,  s1, 
                               h20, s20, 
                               h10, s10, 
                               eta1, eta2, 
                               w=w)
            if not restart:
                stm_calc.initialize(energies, bias = bias)
            
            self.stm_calc = stm_calc
            self.transport_uptodate = True            
            print >> self.log, '#%d:%02d:%02d'\
                     % (T[3], T[4], T[5])+' Done'
            self.log.flush()

    def set_tip_position(self, position):   
        """Positions origin of the tip unit cell at grdpt_c
            in the extended transverse surface unit cell."""         
        p = self.input_parameters
        tip_atom_index = p['tip_atom_index']
        h_c = self.srf.gd.h_c        
        if len(position) == 2:
            position = np.resize(position,3)
            position[2] = 0

        tip_pos_av = self.tip_cell.atoms.get_positions() / Bohr
        tip_pos_av_grpt = self.tip_cell.gd.N_c\
                          * self.tip_cell.atoms.get_scaled_positions()
        srf_pos_av = self.srf.atoms.get_positions() / Bohr
        tip_zmin = tip_pos_av[tip_atom_index, 2]
        srf_zmax = srf_pos_av[:, 2].max()
        extension_c = np.resize(self.srf_cell.ext1,3)
        extension_c[-1] =0
        self.extension_c = extension_c
        #corner of the tip unit cell in the extended grid        
        cell_corner_c = position + self.extension_c\
                      - tip_pos_av_grpt[tip_atom_index]
        cell_corner_c[2]  = (srf_zmax + self.dmin - tip_zmin) / h_c[2]
        cell_corner_c = np.round(cell_corner_c).astype(int)        

        self.tip_position = cell_corner_c + \
                            tip_pos_av_grpt[tip_atom_index] - self.extension_c        
        self.dmin = self.tip_position[2] * h_c[2] - srf_zmax
        self.tip_cell.set_position(cell_corner_c)        

        # Add the tip potential at the respective place in the extended grid
        # of the surface
        size_c = self.tip_cell.gd.n_c
        #current_Vt = self.extd_vt_G.copy()#XXX
        
        current_Vt = self.srf_cell.vt_G.copy()
        current_Vt[cell_corner_c[0]+1:cell_corner_c[0] + 1 + size_c[0],
                   cell_corner_c[1]+1:cell_corner_c[1] + 1 + size_c[1],
                   cell_corner_c[2]+1:cell_corner_c[2] + 1 + size_c[2]] +=\
                    self.tip_cell.vt_G # +1 since cell grid starts at (1,1,1), pbc = 0
        self.current_v=current_Vt
        self.tip_cell_corner_c = cell_corner_c

    def get_V(self, position_c):
        """Returns the overlap hamiltonian for a position of the tip_atom """
        if not self.initialized:
            self.initialize()
        f_iGs = self.srf_cell.f_iGs
        self.set_tip_position(position_c)
        nj = self.nj
        ni = self.ni
        V_ij = np.zeros((nj, ni))
        S_ij = np.zeros((nj, ni))
        vt_G = self.current_v 
        for s in self.srf_cell.functions:
            j1 = s.index
            s.f_iG = f_iGs[j1]
            j2 = j1 + len(s)
            for t, t_kin in zip(self.tip_cell.functions,\
                                self.tip_cell.functions_kin):
                i1 = t.index
                i2 = i1 + len(t)
                V = (s | vt_G | t) 
                S = (s | t )
                if V is None:
                    V = 0
                kin = (s | t_kin)
                if kin is None:
                    kin = 0
                V_ij[j1:j2, i1:i2] += V + kin
                if S is not None:
                    S_ij[j1:j2,i1:i2] += S
            s.f_iG = None
        return V_ij * Hartree - S_ij * self.ediff.real

    def get_transmission(self, position_c):
        V_ts = self.get_V(position_c)       
        T_stm = self.stm_calc.get_transmission(V_ts)
        return T_stm

    def get_current(self, position_c, bias=None):
        if bias == None:
            bias = self.stm_calc.bias
        position_c = tuple(position_c)+(0,)
        V_ts = self.get_V(position_c)
        Is = self.stm_calc.get_current(bias,V_ts)    
        return Is*77466.1509   #units: nA
    
    def get_s(self, position_c):
        self.set_tip_position(position_c)
        S_ij = np.zeros((self.nj, self.ni))
        for s in self.srf_functions:
            j1 = s.index
            j2 = j1 + len(s)
            for t in self.tip_functions:
                i1 = t.index
                i2 = i1 + len(t)
                overlap = (s | t) 
                if overlap is not None:
                    S_ij[j1:j2, i1:i2] += overlap
  
        return S_ij

    def reset(self):
        self.scans = {}

    def scan(self):
        n_c = self.srf.gd.n_c
        scan = np.zeros(tuple(n_c[:2]))
        for x in range(n_c[0]):
             for y in range(n_c[1]):
                 I  = self.get_current([x,y])                       
                 scan[x,y] = I
                 if self.input_parameters['verbose']:
                    print self.tip_cell.position, I
        self.scans['fullscan'] = scan

    def linescan(self, startstop=None):
        if startstop == None:
            start = np.array([0,0])
            stop = self.srf.gd.n_c[:2]-1
        else:
            start = np.asarray(startstop[0])
            stop = np.asarray(startstop[1])
        assert ((self.srf.gd.n_c[:2]-stop)>=0).all()
        v = (stop-start)/np.linalg.norm(stop-start)
        self.v = v
        n_c = self.srf.gd.n_c[:2]
        h_c = self.srf.gd.h_c[:2]
        h = np.linalg.norm(v*h_c)
        n = np.floor(np.linalg.norm((stop-start)*h_c)/h).astype(int) + 1
        linescan_n = np.zeros((n,))
        line = np.arange(n)*Bohr*h
        for i in range(n):
            grpt = start + v*i
            if np.round((grpt%1),5).any(): # Interpolate if nessesary
                C = np.empty((2,2,2)) # find four nearest neighbours
                C[0,0] = np.floor(grpt)
                C[1,0] = C[0,0] + np.array([1,0])
                C[0,1] = C[0,0] + np.array([0,1])
                C[1,1] = C[0,0] + np.array([1,1])
                xd = (grpt%1)[0]
                yd = (grpt%1)[1]
                if not self.scans.has_key('fullscan'):
                    I1 = self.get_current(C[0,0])*(1-xd) \
                        + self.get_current(C[1,0])*xd
                    I2 = self.get_current(C[0,1])*(1-xd) \
                        + self.get_current(C[1,1])*xd
                    I = I1 * (1-yd) + I2*yd
                else:
                    fullscan = self.scans['fullscan']
                    I1 = fullscan[tuple(C[0,0])]*(1-xd) \
                       + fullscan[tuple(C[1,0])]*xd
                    I2 = fullscan[tuple(C[0,1])]*(1-xd) \
                       + fullscan[tuple(C[1,1])]*xd
                    I = I1 * (1-yd) + I2*yd

            else:
                if not self.scans.has_key('fullscan'):
                    I = self.get_current(grpt.astype(int))
                else:                
                    I = self.scans['fullscan'][tuple(grpt.astype(int))]
            linescan_n[i] = I
        self.scans['linescan']=([start,stop],line,linescan_n) 

    def write(self, filename):
        stmc = self.stm_calc
        fd = open(filename,'wb')
        pickle.dump({'p': self.input_parameters,
                     'egft12_emm': [stmc.energies, stmc.gft1_emm, stmc.gft2_emm]
                    },fd,2)
        fd.close()

    def get_dmin(self):
        return self.dmin * Bohr

    def restart(self, filename):
        restart = pickle.load(open(filename))
        p = restart['p']
        self.set(**p)
        print >> self.log, '#Restarting from restart file'
        print >> self.log, '#   bias = ' + str(p['bias'])
        print >> self.log, '#   de = ' + str(p['de'])
        print >> self.log, '#   w = ' + str(p['w'])
        self.initialize_transport(restart = True) 
        self.stm_calc.energies = restart['egft12_emm'][0]
        self.stm_calc.gft1_emm = restart['egft12_emm'][1]    
        self.stm_calc.gft2_emm = restart['egft12_emm'][2]    
        self.stm_calc.bias = p['bias']
        class Dummy:
            def __init__(self, bias):
                self.bias = bias
        self.stm_calc.selfenergy1 = Dummy(p['bias'] * p['w'])    
        self.stm_calc.selfenergy2 = Dummy(p['bias'] * (p['w'] - 1))    
        self.log.flush()

    def hs_from_paw(self): # Do not even try this for larger systems
        p=self.input_parameters
        self.set(**{'hs1':dump_hs(self.tip, 'hs1', 1, 1),
                    'hs2':dump_hs(self.srf, 'hs2', 1, 1),
                    'hs10':dump_lead_hs(self.lead1, 'hs10', 1, 1),
                    'hs20':dump_lead_hs(self.lead2, 'hs20', 1, 1)})

    def mem_estimate(self):
        # basis functions
        M_bfs = 0
        for ft, ft_kin in zip(self.tip_cell.functions, self.tip_cell.functions_kin):
            M_bfs += len(np.ravel(ft.f_iG))
            M_bfs += len(np.ravel(ft_kin.f_iG))
        for f_iG in self.srf_cell.f_iGs.values():
                M_bfs += len(np.ravel(f_iG))
        M_bfs *= 8 / 1024.**2
        print 'LFC-objects:', np.round(M_bfs, 2), 'Mb'
        
        # tip green functions
        es = len(self.stm_calc.energies)
        M_tgft = len(self.stm_calc.h2)**2 * es
        M_tgft *= 2 * 8 / 1024.**2 # complex
        print 'tip green functions:', np.round(M_tgft, 2), 'Mb'

        # surface green function
        M_sgft = len(self.stm_calc.h1)**2 * es
        M_sgft *= 2 * 8 /1024.**2 # complex
        print 'surface green functions:', np.round(M_sgft, 2), 'Mb'
        print 'total:', np.round(M_bfs + M_tgft + M_sgft, 2), 'Mb'

    def plot(self, repeat=[1,1], vmin=None, gd = None):
        import matplotlib
        import pylab
        from pylab import ogrid, imshow, cm, colorbar
        
        repeat = np.asarray(repeat)
        if gd == None:
            gd = self.srf.gd
        
        if self.scans.has_key('fullscan'):
            scan0_iG = self.scans['fullscan']
            shape0 = np.asarray(scan0_iG.shape)    
            scan1_iG = np.zeros(shape0 * repeat)
            
            for i in range(repeat[0]):
                for j in range(repeat[1]):
                    start = np.array([i,j]) * shape0
                    stop = start + shape0 
                    scan1_iG[start[0]:stop[0], start[1]:stop[1]] = scan0_iG

            scan0_iG = scan1_iG
            shape = scan0_iG.shape
            scan_iG = np.zeros(tuple(np.asarray(shape)+1))
            scan_iG[:shape[0],:shape[1]] = scan0_iG
            scan_iG[-1,:shape[1]] = scan0_iG[0,:]
            scan_iG[:,-1] = scan_iG[:,0]

            h = 0.2         
            N_c = np.floor((gd.cell_cv[0,:2] * repeat[0]\
                + gd.cell_cv[1, :2] * repeat[1]) / h).astype(int) 
            ortho_cell_c = np.array(N_c * gd.h_c[:2])
            plot = np.zeros(tuple(N_c))

            # Basis change matrix
            # e -> usual basis {(1,0),(0,1)}
            # o -> basis descrining original original cell
            # n -> basis describing the new cell
            eMo = (gd.cell_cv.T / gd.cell_c * gd.h_c)[:2,:2]
            eMn = np.eye(2) * h           
            oMn = np.dot(np.linalg.inv(eMo), eMn)
            
            for i in range(N_c[0]):
                for j in range(N_c[1]):
                    grpt = np.dot(oMn, [i,j])
                    if (grpt<0).any() or (np.ceil(grpt)>shape).any():
                        plot[i,j] = plot.min() - 1000
                    else: # interpolate
                        C00 = np.floor(grpt).astype(int)
                        C01 = C00.copy()
                        C01[0] += 1
                        C10 = C00.copy()
                        C10[1] += 1
                        C11 = C10.copy()
                        C11[0] += 1
                        
                        x0 = grpt[0] - C00[0]
                        y0 = grpt[1] - C00[1]
                        P0 = scan_iG[tuple(C00)] * (1 - x0)\
                           + scan_iG[tuple(C01)] * x0
                        P1 = scan_iG[tuple(C10)] * (1 - x0)\
                           + scan_iG[tuple(C11)] * x0
                        plot[i,j] = P0 * (1 - y0) + P1 * y0

            plot = plot.T # origin to the lower left corner
            self.scans['interpolated_plot'] = plot
            if vmin == None:
                vmin = scan0_iG.min()
            norm = matplotlib.colors.normalize(vmin=vmin, vmax=scan0_iG.max())
            
            f0 = pylab.figure()
            p0 = f0.add_subplot(111)
            x,y = ogrid[0:plot.shape[0]:1, 0:plot.shape[1]:1]
            extent=[0, plot.shape[1] * h * Bohr,
                    0, plot.shape[0] * h * Bohr]        
            
            p0.set_ylabel('\xc5')
            p0.set_xlabel('\xc5')
            imshow(plot,
                   norm=norm,
                   interpolation='bicubic',
                   origin='lower',
                   cmap=cm.hot,
                   extent=extent)
            cb = colorbar()            
            cb.set_label('I[nA]')
                       
        if self.scans.has_key('linescan'):
            startstop, line, linescan_n = self.scans['linescan']
            start = startstop[0]
            stop = startstop[1]
            f1 = pylab.figure()
            p1 = f1.add_subplot(111)
            p1.plot(line, linescan_n)
            eMo = (gd.cell_cv.T / gd.cell_c)[:2,:2]
            start = np.dot(eMo, start * gd.h_c[:2]) * Bohr
            stop = np.dot(eMo, stop * gd.h_c[:2]) * Bohr
             
            if self.scans.has_key('fullscan'): #Add a line
                p0.plot([start[0], stop[0]], [start[1], stop[1]],'-b')
                p0.set_xlim(tuple(extent[:2]))
                p0.set_ylim(tuple(extent[-2:]))
        pylab.show()

class TipCell:
    def __init__(self, tip, srf):
        self.tip = tip
        self.srf = srf
        self.gd = None
        self.vt_G = None
        self.tip_atom_index = None
        self.functions = []
        self.functions_kin = []

    def initialize(self, tip_indices, tip_atom_index, debug=False):
        self.tip_indices = tip_indices
        self.tip_atom_index =  tip_atom_index
        assert tip_atom_index in tip_indices

        tgd = self.tip.gd
        sgd = self.srf.gd
        tip_atoms = self.tip.atoms.copy()[tip_indices]      
        tip_atoms.pbc = 0
        tip_pos_av = tip_atoms.get_positions().copy()/Bohr
        tip_cell_cv = tgd.cell_cv
        srf_cell_cv = sgd.cell_cv

        tip_zmin = tip_pos_av[tip_atom_index, 2]
        tip_zmin_a = np.zeros(len(tip_indices))

        # size of the simulation cell in the z-direction
        m = 0
        for a, setup in enumerate(self.tip.wfs.setups):
            if a in tip_indices:
                rcutmax = max([phit.get_cutoff() for phit in setup.phit_j])
                tip_zmin_a[m] = tip_pos_av[a, 2] - rcutmax - tip_zmin
                m+=1
        p=2
        zmax_index = np.where(tip_pos_av[:, 2] == tip_pos_av[:, 2].max())[0][0]
        cell_zmin = tip_zmin + tip_zmin_a.min()
        cell_zmax = 2 * tip_pos_av[zmax_index, 2]\
                  - tip_zmin - tip_zmin_a[zmax_index]
        
        if cell_zmax > tgd.cell_c[2] - tgd.h_c[2]:
            cell_zmax = tgd.cell_c[2] - tgd.h_c[2]

        cell_zmin_grpt = np.floor(cell_zmin / tgd.h_c[2] - p).astype(int)
        cell_zmax_grpt = np.floor(cell_zmax / tgd.h_c[2]).astype(int)
        new_sizez = cell_zmax_grpt - cell_zmin_grpt
        self.cell_zmax_grpt = cell_zmax_grpt
        self.cell_zmin_grpt = cell_zmin_grpt

        # If tip and surface cells differ in the xy-plane, 
        # determine the 2d-cell with the smallest area, having lattice vectors 
        # along those vectors describing the 2d-cell belonging to the surface. 

        srf_basis = srf_cell_cv.T / sgd.cell_c
        tip_basis = tip_cell_cv.T / tgd.cell_c

        if (srf_basis - tip_basis).any(): # different unit cells   
            dointerpolate = True
            steps = 500 # XXX crap
            thetas = np.arange(0, pi, pi/steps) 
            areas = np.zeros_like(thetas)
            for i, theta in enumerate(thetas):
                cell = smallestbox(tip_cell_cv, srf_cell_cv, theta)[0]
                area = np.cross(cell[0, :2], cell[1, :2])
                areas[i] = abs(area)
            area_min_index = np.where(areas == areas.min())[0].min()
            theta_min = thetas[area_min_index]
            newcell_cv, origo_c = smallestbox(tip_cell_cv, srf_cell_cv, theta_min)
            tip_pos_av = np.dot(rotate(theta_min), tip_pos_av.T).T + origo_c
            newcell_c = np.array([la.norm(cell_cv[x]) for x in range(3)])
            newsize_c = np.around(newcell_c / sgd.h_c).astype(int)
        elif (sgd.h_c - tgd.h_c).any(): # different grid spacings
            dointerpolate = True
            newcell_cv = tip_cell_cv
            newcell_c = tgd.cell_c
            newsize_c = np.around(newcell_c / sgd.h_c).astype(int)
            theta_min = 0.0
            origo_c = np.array([0,0,0])
        else:
            dointerpolate = False
            newsize_c = tgd.N_c
            vt_G = self.tip.hamiltonian.vt_sG[0]
            vt_G = vt_G[:, :, cell_zmin_grpt:cell_zmax_grpt]
            theta_min = 0.0
            origo_c = np.array([0,0,0])
            self.vt_G = vt_G

        tip_pos_av[:,2] -= cell_zmin_grpt * tgd.h_c[2]
        newsize_c[2] = new_sizez
        newcell_c = (newsize_c + 1) * sgd.h_c 
        newcell_cv = srf_basis * newcell_c
        newgd = GridDescriptor(N_c=newsize_c+1,
                               cell_cv=newcell_cv,
                               pbc_c=False,
                               comm=mpi.serial_comm)
        new_basis = newgd.cell_cv.T / newgd.cell_c

        origo_c += np.dot(new_basis, newgd.h_c)
        tip_pos_av += np.dot(new_basis, newgd.h_c)
        tip_atoms.set_positions(tip_pos_av * Bohr)
        tip_atoms.set_cell(newcell_cv * Bohr)
        self.atoms = tip_atoms
        self.gd = newgd
        
        # quick check
        assert not (np.around(new_basis - srf_basis, 5)).all() 
        assert not (np.around(newgd.h_c - sgd.h_c, 5)).all()

        # add functions
        functions = []
        i=0
        for k, a in enumerate(tip_indices):
            setup = self.tip.wfs.setups[a]
            spos_c = self.atoms.get_scaled_positions()[i]
            if a == tip_atom_index:
                for phit in setup.phit_j:
                    f = AtomCenteredFunctions(self.gd, [phit], spos_c, i)
                    functions.append(f)
                    i += len(f.f_iG)
        self.ni = i
           
        # Apply kinetic energy:
        functions_kin = []
        for f in functions:
            functions_kin.append(f.apply_t())
        
        for f, f_kin in zip(functions, functions_kin):
            f.restrict()
            f_kin.restrict()
 
        self.attach(functions,functions_kin)
  
        if dointerpolate:
            self.interpolate_vt_G(theta_min, origo_c)
    
    def attach(self, functions, functions_kin):
        self.functions = functions
        self.functions_kin = functions_kin
        p0 = {}
        p0_kin = {}        

        for f, f_kin in zip(self.functions, self.functions_kin):       
            p0[f]=f.corner_c.copy() 
            p0_kin[f_kin]=f_kin.corner_c.copy()

        self.p0 = p0
        self.p0_kin = p0_kin
        
    def set_position(self, position_c):
        self.position = position_c
        for f, f_kin in zip(self.functions, self.functions_kin):
            f.corner_c =  position_c + self.p0[f]
            f_kin.corner_c = position_c + self.p0_kin[f_kin]
    
    def shift_potential(self, shift):
        pass

    def interpolate_vt_G(self, theta_min, origo_c):
        """Interpolates the effective potential of the tip calculation onto
           the grid of the simulation cell for the tip.       

           The transformation iMj maps a point from grid 'j' to
           a point on grid 'i',  j --> i.
        
           Definitions:
           e - 'natural' grid, {(1, 0, 0), (0, 1, 0), (0, 0, 1)}       
           o - grid of the original tip calculation
           r - rotated grid of the original tip calculation
           n - grid of the tip simulation cell.

           Outside the unitcell of the original tip calculation 
           the effective potential is set to zero"""

        vt_G0 = self.tip.hamiltonian.vt_sG[0]
        vt_G0 = vt_G0[:, :, self.cell_zmin_grpt:self.cell_zmax_grpt]
        tgd = self.tip.gd
        newgd = self.gd
        shape0 = vt_G0.shape
        tip_basis = tgd.cell_cv.T / tgd.cell_c
        new_basis = newgd.cell_cv.T / newgd.cell_c

        eMo = tip_basis * tgd.h_c
        eMr = np.dot(rotate(theta_min), eMo)
        eMn = new_basis * newgd.h_c
        rMn = np.dot(la.inv(eMr), eMn)

        vt_G = newgd.zeros()
        shape = vt_G.shape
        self.shape2 = shape
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    gpt_n = [i, j, k]
                    gpt_r = np.dot(rMn, gpt_n) - np.dot(la.inv(eMr), origo_c)
                    if (gpt_r < 0).any() or (np.ceil(gpt_r) > tgd.n_c).any():
                        vt_G[i,j,k] =  0
                    else: # trilinear interpolation
                        C000 = np.floor(gpt_r).astype(int)
                        z00 = gpt_r[2] - C000[2]
                        C001 = C000.copy()
                        C001[2] += 1
                        C100 = C000.copy()
                        C100[0] += 1
                        C101 = C100.copy()
                        C101[2] += 1
                        C010 = C000.copy()
                        C010[1] += 1
                        C011 = C010.copy()
                        C011[2] += 1
                        C110 = C000.copy()
                        C110[:2] += 1
                        C111 = C110.copy()
                        C111[2] += 1
                        x0 = gpt_r[0] - C000[0]
                        y0 = gpt_r[1] - C000[1]
                        C = np.zeros((4,2))
                        C1 = np.array([[vt_G0[tuple(C000)], vt_G0[tuple(C001 % shape0)]],
                                      [ vt_G0[tuple(C010 % shape0)], vt_G0[tuple(C011 % shape0)]]])
                        C2 = np.array([[vt_G0[tuple(C100 % shape0)], vt_G0[tuple(C101 % shape0)]],
                                      [ vt_G0[tuple(C110 % shape0)], vt_G0[tuple(C111 % shape0)]]])
                        Z = np.array([1 - z00, z00])
                        X = np.array([1 - x0, x0])
                        Y = np.array([1 - y0, y0])
                        Q = np.zeros((2, 2))
                        Q[:,0]=np.dot(C1, Z)
                        Q[:,1]=np.dot(C2, Z)
                        F2 = dots(Y, Q, X)
                        vt_G[i, j, k] = F2
        self.vt_G = vt_G

class SrfCell:
    def __init__(self, srf):
        self.srf = srf
        self.functions = []

    def initialize(self, tip_cell, srf_indices):
        self.srf_indices = srf_indices
        # determine the extended unitcell
        srf_vt_G = self.srf.hamiltonian.vt_sG[0]
        tip = tip_cell
        tip_atom_index = tip.tip_atom_index
        spos_ac = tip.atoms.get_scaled_positions()
        tip_atom_spos = spos_ac[tip_atom_index][:2]
        tgd = tip.gd
        sgd = self.srf.gd
        tip_cell_cv = tgd.cell_cv[:2, :2]
        tip_cell_c = tgd.cell_c[:2]
        tip_basis = tip_cell_cv.T / tip_cell_c
        srf_cell_cv = sgd.cell_cv[:2, :2]
        srf_cell_c = sgd.cell_c[:2]
        srf_basis = tip_cell_cv.T / tip_cell_c
        assert not (np.round(tgd.h_c - sgd.h_c, 5)).all()
        assert not (np.round(tip_basis - srf_basis, 5)).all()
        extension1_c = tip_atom_spos * tip_cell_c / srf_cell_c
        extension2_c = (1 - tip_atom_spos) * tip_cell_c / srf_cell_c
        ext1_c = np.ceil(extension1_c * sgd.n_c[:2]).astype(int)
        ext2_c = np.ceil(extension2_c * sgd.n_c[:2]).astype(int)

        srf_shape  = sgd.n_c[:2]
        extension1 = ext1_c / srf_shape.astype(float)
        extension2 = ext2_c / srf_shape.astype(float)
        newsize_c = ext1_c + ext2_c + sgd.n_c[:2]
        vt_G = np.zeros(tuple(newsize_c) + (srf_vt_G.shape[2],))
 
        intexa = ext1_c / srf_shape[:2]
        rest1 = ext1_c % srf_shape[:2]
        intexb = ext2_c / srf_shape[:2]
        rest2 = ext2_c % srf_shape[:2]
        print intexa, intexb, rest1, rest2

        for n in range(intexa[0]+intexb[0] + 1 ):
            for m in range(intexa[1] + intexb[1] + 1):
                vt_G[rest1[0] + n * srf_shape[0]:\
                     rest1[0] + (n + 1) * srf_shape[0],\
                     rest1[1] + m * srf_shape[1]:\
                     rest1[1] + (m + 1) * srf_shape[1]] = srf_vt_G

        if rest2[1] == 0:
            rest2[1] += 1
        if rest2[0] == 0:
            rest2[0] += 1

        vt_G[:rest1[0], rest1[1]: -rest2[1]]\
             = vt_G[-rest1[0] - rest2[0]:-rest2[0], rest1[1]:-rest2[1]]
        vt_G[-rest2[0]:, rest1[1]:-rest2[1]]\
             = vt_G[rest1[0]:rest1[0] + rest2[0], rest1[1]:-rest2[1]]
        vt_G[:, :rest1[1]] = vt_G[:, -rest2[1] - rest1[1]:-rest2[1]]
        vt_G[:, -rest2[1]:] = vt_G[:, rest1[1]:rest1[1]+rest2[1]]
        
        self.vt_G = vt_G
        newsize_c = np.resize(newsize_c, 3)
        newsize_c[2] = sgd.n_c[2]
        newcell_cv = (newsize_c+1) * sgd.cell_cv.T / sgd.cell_c * sgd.h_c
        newgd = GridDescriptor(N_c=newsize_c + 1,
                               cell_cv=newcell_cv,
                               pbc_c=False,
                               comm=mpi.serial_comm)

        self.gd = newgd
        srf_atoms = self.srf.atoms.copy()[srf_indices]
        self.atoms = srf_atoms
        
        # add functions
        j = 0
        for a in srf_indices:
            setup = self.srf.wfs.setups[a]
            spos_c = self.srf.atoms.get_scaled_positions()[a]
            for phit in setup.phit_j:
                f = AtomCenteredFunctions(self.srf.gd, [phit], spos_c, j)
                self.functions.append(f)
                j += len(f.f_iG)
        self.nj = j
        # shift corners so that the origin now is the one of the extended surface cell
        for f in self.functions:
            f.corner_c[:2] += ext1_c #XXX maybe not + 1 ???
        self.ext1 = ext1_c

        # Add an appropriate number of periodic images.
        # Translation vectors:
        Rs =np.array([[0,  1],
                     [1,   1],
                     [1,   0],
                     [1,  -1],
                     [0,  -1],
                     [-1, -1],
                     [-1,  0],
                     [-1,  1]])

        origo = np.array([0, 0])
        list = []
        f_iGs = {}
        for f in self.functions:
            f_iGs[f.index] = f.f_iG
            f.f_iG = None
            list.append(f)
            for R in Rs:
                n = 0
                add_function = True
                while add_function == True:
                    n += 1
                    newcorner_c = f.corner_c[:2] +  n * R * sgd.n_c[:2]
                    start_c = np.maximum(newcorner_c, origo)
                    stop_c = np.minimum(newcorner_c + f.size_c[:2],
                                        newgd.n_c[:2])
                    if (start_c < stop_c).all():
                        newcorner_c = np.resize(newcorner_c, 3)
                        newcorner_c[2] = f.corner_c[2]
                        newf = LocalizedFunctions(f.gd, f_iGs[f.index],
                                                  corner_c=newcorner_c,
                                                  index=f.index,
                                                  vt_G=f.vt_G)
                        newf.f_iG = None
                        newf.set_phase_vector(n * R)
                        list.append(newf)
                    else:
                        add_function = False
        self.functions = list
        self.f_iGs = f_iGs
        self.atoms = srf_atoms

    def shift_potential(self, shift):
        pass

def dump_hs(calc, filename, return_hs=False, restart = False):
    """Pickle LCAO - Hamiltonian and overlap matrix for a tip or surface
    calculation.
    
    use restart=True if restarting from a gpw file. The energy scale is set
    to zero at the fermi level and periodic boundary conditions are removed 
    in the z-direction.
    """
    atoms = calc.get_atoms()
    atoms.set_calculator(calc)

    if restart:
       calc.initialize_positions(atoms)

    efermi = calc.get_fermi_level()

    if w.rank == 0:
        calc.wfs.eigensolver.calculate_hamiltonian_matrix(
                 calc.hamiltonian, calc.wfs, calc.wfs.kpt_u[0]) # for now gamma only
        s_mm = calc.wfs.S_qMM[0].copy() 
        h_mm = calc.wfs.eigensolver.H_MM * Hartree
        h_mm -= s_mm * efermi
        tri2full(s_mm)
        tri2full(h_mm)
        remove_pbc(atoms, h_mm, s_mm, 2)
        fd = open(filename + '.pckl', 'wb')        
        pickle.dump((h_mm, s_mm), fd)
        fd.close()
        if return_hs:
            return (h_mm, s_mm)


def dump_lead_hs(calc, filename, return_hs=False, restart=False):
    """Pickle real space LCAO - Hamiltonian and overlap matrix for a 
    periodic lead calculation.

    use restart=True if restarting from a gpw file.
    The energy scale is set to zero at fermi level.
    """
    atoms = calc.get_atoms()
    atoms.set_calculator(calc)
 
    if restart:
       calc.initialize_positions(atoms)

    efermi = calc.get_fermi_level()
    h_smm, s_mm = get_lead_lcao_hamiltonian(calc, direction='z')
    h_mm = h_smm[0] - efermi * s_mm
    fd = open(filename + '.pckl', 'wb')
    pickle.dump((h_mm, s_mm), fd, 2)
    fd.close()
    if return_hs:
        return (h_mm, s_mm)
 
def intersection(l1, l2):
    """Intersection (x, y, t) between two lines.
    
    Two points on each line have to be specified.
    """
    a1 = l1[0]
    b1 = l1[1]

    a2 = l2[0]
    b2 = l2[1]

    A = np.zeros((2,2))
    A[:,0] = b1-a1
    A[:,1] = a2-b2

    if np.round(np.linalg.det(A),5) == 0: #parallel lines
        return None
    
    r = a2 - a1
    t = np.dot(la.inv(A), r.T)
    xy = a2 + t[1] * (b2 - a2)
    return (list(xy), t)

def rotate(theta):
    return np.array([[cos(theta), sin(theta), 0],
                    [-sin(theta), cos(theta), 0],
                    [0, 0, 1]])

def unravel2d(data, shape):
    pass

def smallestbox(cell1, cell2, theta, plot=False):
    """Determines the smallest 2d unit cell which encloses cell1 rotated at 
    an angle theta around the z-axis and which has lattice vectors parallel 
    to those of cell2."""

    ct = cell1[:2,:2] * Bohr
    cs = cell2[:2,:2] * Bohr
    v3 = [cs[0],cs[1]]
    lsrf = [np.array([0.,0.]),v3[1],v3[0]+v3[1],v3[0],np.array([0,0])]
    cs = cs/np.array([la.norm(cs[0]),la.norm(cs[1])])
    v4 = [ct[0],ct[1]]
    lct = [np.array([0.,0.]),v4[1],v4[0]+v4[1],v4[0],np.array([0,0])]
    new_ct = np.dot(R(theta), ct.T).T
    v1 = [new_ct[0],new_ct[1]]
    v2 = [cs[0],cs[1]]
    l1 = [np.array([0.,0.]),v1[1],v1[0]+v1[1],v1[0],np.array([0,0])]

    sides1 = []
    for j in range(4):
        intersections = 0
        for i in range(4):
            line = (l1[j], l1[j]+v2[0])
            side = (l1[i], l1[i+1])  
            test = intersection(line,side)
            if test is not None:
                t = test[1]
                if np.round(t[1],5)<=1 and np.round(t[1],5)>=0:
                    intersections += 1
        if intersections == 2:
            sides1.append(line)

    sides2 = []
    for j in range(4):
        intersections = 0
        for i in range(4):
            line = (l1[j],l1[j]+v2[1])
            side = (l1[i], l1[i+1])  
            test = intersection(line,side)
            if test is not None:
                t = test[1]
                if np.round(t[1],5)<=1 and np.round(t[1],5)>=0:
                    intersections += 1
        if intersections == 2:
            sides2.append(line)
    corners = []
    for side1 in sides1:
        for side2 in sides2:
            corner=np.round(intersection(side1,side2),5)
            if corner is not None:
                corners.append(list(corner[0]))
    for i in range(len(corners)):
        for j in range(2):
            if np.round(corners[i][j],5) == 0:
                corners[i][j] = abs(corners[i][j])

    corners.sort()   
    if len(corners)==8:
        corners = [corners[0],corners[2],corners[4],corners[6]]
    if len(corners)==16:
        corners = [corners[0],corners[4],corners[8],corners[12]]
    origo = np.array(corners.pop(0))
    end =np.array(corners.pop(-1))
    pa = corners[0]
    pb = corners[1]
    if pa[1]/pa[0] < pb[1]/pb[0]:
        lat1 = pa-origo
        lat2 = pb-origo
    else:
        lat1 = pb-origo
        lat2 = pa-origo
    l1-=origo
    area = np.cross(lat1,lat2)

    '''
    if plot:
        import pylab
        f0 = pylab.figure()
        p0 = f0.add_subplot(111)
        p0.set_title(str(theta))
        ltest = [np.array([0,0]),lat1,lat2+lat1,lat2,np.array([0,0])]
        for i in range(4):
            start = l1[i]
            stop = l1[i+1]
            p0.plot([start[0],stop[0]],[start[1],stop[1]],'-g')
        for i in range(4):
            start = lsrf[i]
            stop = lsrf[i+1]
            p0.plot([start[0],stop[0]],[start[1],stop[1]],'-b')
        for i in range(4):
            start = lct[i]
            stop = lct[i+1]
            p0.plot([start[0],stop[0]],[start[1],stop[1]],'--y')
        for i in range(4):
            start = ltest[i]
            stop = ltest[i+1]
            p0.plot([start[0],stop[0]],[start[1],stop[1]],'-r')
        p0.set_xlim([-1,max(lsrf[2][0],(end-origo)[0])+1])
        p0.set_ylim([-1,max(lsrf[2][1], (end-origo)[1])+1])
        pylab.show()
    cell = np.zeros((3,3))
    cell[2,2] = 1
    cell[0,:2] = lat1
    cell[1,:2] = lat2
    origo = np.zeros(3,)
    origo[:2] = l1[0]
    return cell/Bohr, origo/Bohr
    '''

