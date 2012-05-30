import numpy as np
from math import sqrt, pi
import pickle
from ase.units import Hartree, Bohr
from gpaw.mpi import rank
from gpaw.response.chi import CHI


class DF(CHI):
    """This class defines dielectric function related physical quantities."""

    def __init__(self,
                 calc=None,
                 nbands=None,
                 w=None,
                 q=None,
                 eshift=None,
                 ecut=10.,
                 smooth_cut=None,
                 density_cut=None,
                 G_plus_q=False,
                 eta=0.2,
                 rpad=np.array([1, 1, 1]),
                 vcut=None,
                 ftol=1e-7,
                 txt=None,
                 xc='ALDA',
                 print_xc_scf=False,
                 hilbert_trans=True,
                 full_response=False,
                 optical_limit=False,
                 comm=None,
                 kcommsize=None):

        CHI.__init__(self, calc=calc, nbands=nbands, w=w, q=q, eshift=eshift,
                     ecut=ecut, smooth_cut=smooth_cut, density_cut=density_cut,
                     G_plus_q=G_plus_q, eta=eta, rpad=rpad, vcut=vcut,
                     ftol=ftol, txt=txt, xc=xc, hilbert_trans=hilbert_trans,
                     full_response=full_response, optical_limit=optical_limit,
                     comm=comm, kcommsize=kcommsize)

        self.df_flag = False
        self.print_bootstrap = print_xc_scf
        self.df1_w = None  # NLF RPA
        self.df2_w = None  # LF RPA
        self.df3_w = None  # NLF ALDA
        self.df4_w = None  # LF ALDA

    def get_dielectric_matrix(self, xc='RPA', overwritechi0=False):
        
	if self.chi0_wGG is None:
            self.initialize()
            self.calculate()
        else:
            pass  # read from file and re-initializing .... need to be implemented

        if overwritechi0:
            dm_wGG = self.chi0_wGG
        else:
            dm_wGG = np.zeros_like(self.chi0_wGG)
        if xc == 'RPA':
            Kc_G = np.zeros(self.npw)
            for iG in range(self.npw):
                qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
                Kc_G[iG] = 1. / np.sqrt(np.dot(qG, qG))
            Kc_G *= np.sqrt(4 * pi)

            self.printtxt('Use RPA.')
            for iw in range(self.Nw_local):
                for iG in range(self.npw):
                    dmtmp_G = - Kc_G[iG] * Kc_G[:] * self.chi0_wGG[iw,iG,:]
                    dmtmp_G[iG] += 1
                    dm_wGG[iw,iG,:] = dmtmp_G[:]

        elif xc == 'ALDA':
            self.printtxt('Use ALDA kernel.')
            # E_LDA = 1 - v_c chi0 (1-fxc chi0)^-1
            # http://prb.aps.org/pdf/PRB/v33/i10/p7017_1 eq. 4
            tmp_GG = np.eye(self.npw, self.npw)
            A_wGG = self.chi0_wGG.copy()
            for iw in range(self.Nw_local):
                A_wGG[iw] = np.dot(self.chi0_wGG[iw], np.linalg.inv(tmp_GG - np.dot(self.Kxc_sGG[0], self.chi0_wGG[iw])))
    
            for iw in range(self.Nw_local):
                dm_wGG[iw] = tmp_GG - self.Kc_GG * A_wGG[iw]                

        return dm_wGG


    def get_inverse_dielectric_matrix(self,xc='RPA'):

        dm_wGG = self.get_dielectric_matrix(xc=xc)
        dminv_wGG = np.zeros_like(dm_wGG)
        for iw in range(self.Nw_local):
            dminv_wGG[iw] = np.linalg.inv(dm_wGG[iw])
        return dminv_wGG


    def get_chi(self, xc='RPA'):
        """Solve Dyson's equation."""

	if self.chi0_wGG is None:
            self.initialize()
            self.calculate()
        else:
            pass # read from file and re-initializing .... need to be implemented

        kernel_GG = np.zeros((self.npw, self.npw), dtype=complex)
        chi_wGG = np.zeros_like(self.chi0_wGG)

        # Coulomb kernel
        for iG in range(self.npw):
            qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
            kernel_GG[iG,iG] = 4 * pi / np.dot(qG, qG)
            
        if xc == 'ALDA':
            kernel_GG += self.Kxc_sGG[0]

        for iw in range(self.Nw_local):
            tmp_GG = np.eye(self.npw, self.npw) - np.dot(self.chi0_wGG[iw], kernel_GG)
            chi_wGG[iw] = np.dot(np.linalg.inv(tmp_GG) , self.chi0_wGG[iw])

        return chi_wGG
    

    def get_dielectric_function(self, xc='RPA'):
        """Calculate the dielectric function. Returns df1_w and df2_w.

        Parameters:

        df1_w: ndarray
            Dielectric function without local field correction.
        df2_w: ndarray
            Dielectric function with local field correction.
        """

        if self.df_flag is False:
            dm_wGG = self.get_dielectric_matrix(xc=xc)

            Nw_local = dm_wGG.shape[0]
            dfNLF_w = np.zeros(Nw_local, dtype = complex)
            dfLFC_w = np.zeros(Nw_local, dtype = complex)
            df1_w = np.zeros(self.Nw, dtype = complex)
            df2_w = np.zeros(self.Nw, dtype = complex)

            for iw in range(Nw_local):
                tmp_GG = dm_wGG[iw]
                dfLFC_w[iw] = 1. / np.linalg.inv(tmp_GG)[0, 0]
                dfNLF_w[iw] = tmp_GG[0, 0]

            self.wcomm.all_gather(dfNLF_w, df1_w)
            self.wcomm.all_gather(dfLFC_w, df2_w)

            if xc == 'RPA':
                self.df1_w = df1_w
                self.df2_w = df2_w
            elif xc=='ALDA' or xc=='ALDA_X':
                self.df3_w = df1_w
                self.df4_w = df2_w                

        if xc == 'RPA':
            return self.df1_w, self.df2_w
        elif xc == 'ALDA' or xc=='ALDA_X':
            return self.df3_w, self.df4_w


    def get_surface_response_function(self, z0=0., filename='surf_EELS'):
        """Calculate surface response function."""

	if self.chi0_wGG is None:
            self.initialize()
            self.calculate()


        g_w2 = np.zeros((self.Nw,2), dtype=complex)
        assert self.acell_cv[0,2] == 0. and self.acell_cv[1,2] == 0.

        Nz = self.nG[2] # number of points in z direction
        tmp = np.zeros(Nz, dtype=int)
        nGz = 0         # number of G_z 
        for i in range(self.npw):
            if self.Gvec_Gc[i, 0] == 0 and self.Gvec_Gc[i, 1] == 0:
                tmp[nGz] = self.Gvec_Gc[i, 2]
                nGz += 1
        assert (np.abs(self.Gvec_Gc[:nGz, :2]) < 1e-10).all()

        for id, xc in enumerate(['RPA', 'ALDA']):
            chi_wGG = self.get_chi(xc=xc)
    
            # The first nGz are all Gx=0 and Gy=0 component
            chi_wgg_LFC = chi_wGG[:, :nGz, :nGz]
            del chi_wGG
            chi_wzz_LFC = np.zeros((self.Nw_local, Nz, Nz), dtype=complex)        
    
            # Fourier transform of chi_wgg to chi_wzz
            Gz_g = tmp[:nGz] * self.bcell_cv[2,2]
            z_z = np.linspace(0, self.acell_cv[2,2]-self.h_cv[2,2], Nz)
            phase1_zg = np.exp(1j  * np.outer(z_z, Gz_g))
            phase2_gz = np.exp(-1j * np.outer(Gz_g, z_z))
    
            for iw in range(self.Nw_local):
                chi_wzz_LFC[iw] = np.dot(np.dot(phase1_zg, chi_wgg_LFC[iw]), phase2_gz)
            chi_wzz_LFC /= self.acell_cv[2,2]        
    
            # Get surface response function
    
            z_z -= z0 / Bohr
            q_v = np.dot(self.q_c, self.bcell_cv)
            qq = sqrt(np.inner(q_v, q_v))
            phase1_1z = np.array([np.exp(qq*z_z)])
            phase2_z1 = np.exp(qq*z_z)
    
            tmp_w = np.zeros(self.Nw_local, dtype=complex)        
            for iw in range(self.Nw_local):
                tmp_w[iw] = np.dot(np.dot(phase1_1z, chi_wzz_LFC[iw]), phase2_z1)[0]            
    
            tmp_w *= -2 * pi / qq * self.h_cv[2,2]**2        
            g_w = np.zeros(self.Nw, dtype=complex)
            self.wcomm.all_gather(tmp_w, g_w)
            g_w2[:, id] = g_w
    
        if rank == 0:
            f = open(filename,'w')
            for iw in range(self.Nw):
                energy = iw * self.dw * Hartree
                print >> f, energy, np.imag(g_w2[iw, 0]), np.imag(g_w2[iw, 1])
            f.close()

        # Wait for I/O to finish
        self.comm.barrier()


    def check_sum_rule(self, df1_w=None, df2_w=None):
        """Check f-sum rule."""

	if df1_w is None:
            df1_w = self.df1_w
            df2_w = self.df2_w

        N1 = N2 = 0
        for iw in range(self.Nw):
            w = iw * self.dw
            N1 += np.imag(df1_w[iw]) * w
            N2 += np.imag(df2_w[iw]) * w
        N1 *= self.dw * self.vol / (2 * pi**2)
        N2 *= self.dw * self.vol / (2 * pi**2)

        self.printtxt('')
        self.printtxt('Sum rule for ABS:')
        nv = self.nvalence
        self.printtxt('Without local field: N1 = %f, %f  %% error' %(N1, (N1 - nv) / nv * 100) )
        self.printtxt('Include local field: N2 = %f, %f  %% error' %(N2, (N2 - nv) / nv * 100) )

        N1 = N2 = 0
        for iw in range(self.Nw):
            w = iw * self.dw
            N1 -= np.imag(1/df1_w[iw]) * w
            N2 -= np.imag(1/df2_w[iw]) * w
        N1 *= self.dw * self.vol / (2 * pi**2)
        N2 *= self.dw * self.vol / (2 * pi**2)
                
        self.printtxt('')
        self.printtxt('Sum rule for EELS:')
        nv = self.nvalence
        self.printtxt('Without local field: N1 = %f, %f  %% error' %(N1, (N1 - nv) / nv * 100) )
        self.printtxt('Include local field: N2 = %f, %f  %% error' %(N2, (N2 - nv) / nv * 100) )


    def get_macroscopic_dielectric_constant(self):
        """Calculate macroscopic dielectric constant. Returns eM1 and eM2

        Macroscopic dielectric constant is defined as the real part of dielectric function at w=0.
        
        Parameters:

        eM1: float
            Dielectric constant without local field correction. (RPA, ALDA)
        eM2: float
            Dielectric constant with local field correction.

        """

        eM1 = np.zeros(2)
        eM2 = np.zeros(2)
        for id, xc in enumerate(['RPA', 'ALDA']):
            df1, df2 = self.get_dielectric_function(xc=xc)
            eM1[id], eM2[id] = np.real(df1[0]), np.real(df2[0])
        self.printtxt('')
        self.printtxt('Macroscopic dielectric constant:')
        self.printtxt('    Without local field (RPA, ALDA): %f, %f' %(eM1[0], eM1[1]) )
        self.printtxt('    Include local field (RPA, ALDA): %f, %f' %(eM2[0], eM2[1]) )        
            
        return eM1, eM2


    def get_absorption_spectrum(self, filename='Absorption.dat'):
        """Calculate optical absorption spectrum. By default, generate a file 'Absorption.dat'.

        Optical absorption spectrum is obtained from the imaginary part of dielectric function.
        """

        df1, df2 = self.get_dielectric_function(xc='RPA')
        if self.xc == 'ALDA':
            df3, df4 = self.get_dielectric_function(xc='ALDA')
        if self.xc is 'ALDA_X':
            df3, df4 = self.get_dielectric_function(xc='ALDA_X')

        Nw = df1.shape[0]

        if self.xc == 'Bootstrap':
            from gpaw.response.fxc import Bootstrap
            Kc_GG = np.zeros((self.npw, self.npw))
            for iG in range(self.npw):
                qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
                Kc_GG[iG,iG] = 4 * pi / np.dot(qG, qG)

            from gpaw.mpi import world
            assert self.wcomm.size == world.size
            df3 = Bootstrap(self.chi0_wGG, Nw, Kc_GG, self.printtxt, self.print_bootstrap)

        if rank == 0:
            f = open(filename,'w')
            for iw in range(Nw):
                energy = iw * self.dw * Hartree
                if self.xc == 'RPA':
                    print >> f, energy, np.real(df1[iw]), np.imag(df1[iw]), \
                          np.real(df2[iw]), np.imag(df2[iw])
                elif self.xc == 'ALDA':
                    print >> f, energy, np.real(df1[iw]), np.imag(df1[iw]), \
                      np.real(df2[iw]), np.imag(df2[iw]), \
                      np.real(df3[iw]), np.imag(df3[iw]), \
                      np.real(df4[iw]), np.imag(df4[iw])
                elif self.xc == 'Bootstrap':
                    print >> f, energy, np.real(df1[iw]), np.imag(df1[iw]), \
                      np.real(df2[iw]), np.imag(df2[iw]), \
                      np.real(df3[iw]), np.imag(df3[iw])
            f.close()

        # Wait for I/O to finish
        self.comm.barrier()


    def get_EELS_spectrum(self, filename='EELS.dat'):
        """Calculate EELS spectrum. By default, generate a file 'EELS.dat'.

        EELS spectrum is obtained from the imaginary part of the inverse of dielectric function.
        """

        # calculate RPA dielectric function
        df1, df2 = self.get_dielectric_function(xc='RPA')
        if self.xc == 'ALDA':
            df3, df4 = self.get_dielectric_function(xc='ALDA')
        Nw = df1.shape[0]

        if rank == 0:
            f = open(filename,'w')
            for iw in range(self.Nw):
                energy = iw * self.dw * Hartree
                if self.xc == 'RPA':
                    print >> f, energy, -np.imag(1./df1[iw]), -np.imag(1./df2[iw])
                elif self.xc == 'ALDA':
                    print >> f, energy, -np.imag(1./df1[iw]), -np.imag(1./df2[iw]), \
                       -np.imag(1./df3[iw]), -np.imag(1./df4[iw])
            f.close()

        # Wait for I/O to finish
        self.comm.barrier()


    def get_jdos(self, f_skn, e_skn, kd, kq, dw, Nw, sigma):
        """Calculate Joint density of states"""

        JDOS_w = np.zeros(Nw)
        nkpt = kd.nbzkpts
        nbands = f_skn[0].shape[1]

        for k in range(nkpt):
            print k
            ibzkpt1 = kd.bz2ibz_k[k]
            ibzkpt2 = kd.bz2ibz_k[kq[k]]
            for n in range(nbands):
                for m in range(nbands):
                    focc = f_skn[0][ibzkpt1, n] - f_skn[0][ibzkpt2, m]
                    w0 = e_skn[0][ibzkpt2, m] - e_skn[0][ibzkpt1, n]
                    if focc > 0 and w0 >= 0:
                        w0_id = int(w0 / dw)
                        if w0_id + 1 < Nw:
                            alpha = (w0_id + 1 - w0/dw) / dw
                            JDOS_w[w0_id] += focc * alpha
                            alpha = (w0/dw-w0_id) / dw
                            JDOS_w[w0_id+1] += focc * alpha
                            
        w = np.arange(Nw) * dw * Hartree

        return w, JDOS_w


    def calculate_induced_density(self, q, w):
        """ Evaluate induced density for a certain q and w.

        Parameters:

        q: ndarray
            Momentum tranfer at reduced coordinate.
        w: scalar
            Energy (eV).
        """

        if type(w) is int:
            iw = w
            w = self.wlist[iw] / Hartree
        elif type(w) is float:
            w /= Hartree
            iw = int(np.round(w / self.dw))
        else:
            raise ValueError('Frequency not correct !')

        self.printtxt('Calculating Induced density at q, w (iw)')
        self.printtxt('(%f, %f, %f), %f(%d)' %(q[0], q[1], q[2], w*Hartree, iw))

        # delta_G0
        delta_G = np.zeros(self.npw)
        delta_G[0] = 1.

        # coef is (q+G)**2 / 4pi
        coef_G = np.zeros(self.npw)
        for iG in range(self.npw):
            qG = np.dot(q + self.Gvec_Gc[iG], self.bcell_cv)
            coef_G[iG] = np.dot(qG, qG)
        coef_G /= 4 * pi

        # obtain chi_G0(q,w)
        dm_wGG = self.get_RPA_dielectric_matrix()
        tmp_GG = dm_wGG[iw]
        del dm_wGG
        chi_G = (np.linalg.inv(tmp_GG)[:, 0] - delta_G) * coef_G

        gd = self.gd
        r = gd.get_grid_point_coordinates()

        # calculate dn(r,q,w)
        drho_R = gd.zeros(dtype=complex)
        for iG in range(self.npw):
            qG = np.dot(q + self.Gvec_Gc[iG], self.bcell_cv)
            qGr_R = np.inner(qG, r.T).T
            drho_R += chi_G[iG] * np.exp(1j * qGr_R)

        # phase = sum exp(iq.R_i)
        # drho_R /= self.vol * nkpt / phase
        return drho_R


    def get_induced_density_z(self, q, w):
        """Get induced density on z axis (summation over xy-plane). """

        drho_R = self.calculate_induced_density(q, w)

        drho_z = np.zeros(self.nG[2],dtype=complex)
#        dxdy = np.cross(self.h_c[0], self.h_c[1])

        for iz in range(self.nG[2]):
            drho_z[iz] = drho_R[:,:,iz].sum()

        return drho_z

    def get_eigenmodes(self,filename = None, chi0 = None, calc = None, dm = None, xc = 'RPA', sum = None, vcut = None, checkphase = False):
        """
        Calculate the plasmonic eigenmodes by finding eigenvectors of the dielectric matrix.
        If chi0_wGG is given, the RPA dielectric matrix is calculated.  

        Parameters:

        filename:  pckl file
                   output from response calculation.
         
        chi0:  gpw file
               chi0_wGG from response calculation. 

        dm:  gpw file
             dielectric matrix from response calculation  
        
        calc:  gpaw calculator instance
               ground state calculator used in response calculation.
               Wavefunctions only needed if chi0 is calculated from scratch

        sum:  str
              '2D': sum in the x and y directions
              '1D': To be implemented
           
        """
        from gpaw import GPAW
        self.read(filename)      
        self.calc = GPAW(calc,txt=None)
        calc = self.calc
        self.w_w = np.linspace(0, self.dw * (self.Nw - 1)*Hartree, self.Nw)
        self.vcut = vcut
        self.xc = xc
        
        if chi0 is None and dm is None: # calculate chi0 from scratch            
            self.initialize()
            self.calculate()    
            self.w_w = self.w_w * Hartree
            self.dm_wGG = self.get_dielectric_matrix()        
        else: # read from file and reinitialize
            from gpaw.response.parallel import par_read          
            from gpaw.response.kernel import calculate_Kc
            self.Nw_local = self.Nw
            self.nspins = 1
            self.gd = self.calc.wfs.gd
            if dm is not None:
                self.dm_wGG = par_read(dm, 'dm')
            else:
                self.chi0_wGG = par_read(chi0, 'chi0_wGG')
                self.Kc_GG = calculate_Kc(self.q_c, self.Gvec_Gc,
                                          self.acell_cv, self.bcell_cv,
                                          self.calc.atoms.pbc, self.optical_limit,
                                          self.vcut, nonsymmetric = True)
                          
                if self.xc == 'ALDA' or self.xc == 'ALDA_X':
                    # from gpaw.response.kernel import calculate_Kxc
                    self.printtxt('ALDA not implemented')
                    return
                self.dm_wGG = self.get_dielectric_matrix(xc= xc) 
               

        self.pbc = self.calc.atoms.pbc
        q = self.q_c
        gd = self.gd
        r = gd.get_grid_point_coordinates()
        rrr = r*Bohr 
        pbc = self.calc.atoms.pbc
        
        eig_0 = np.array([], dtype = complex)
        eig_left = np.array([], dtype = complex)
        eig_right = np.array([], dtype = complex)
        vec_modes = np.zeros([1, self.npw], dtype = complex)
        vec_modes_dual = np.zeros([1, self.npw], dtype = complex)
        eig_all = np.zeros([1, self.npw], dtype = complex)
        eig_dummy = np.zeros([1, self.npw], dtype = complex)
        v_dummy = np.zeros([1, self.npw], dtype = complex)
        vec_dummy = np.zeros([1, self.npw], dtype = complex)
        w_0 = np.array([]) 
        w_left = np.array([])
        w_right = np.array([])
     
        if sum == '2D':
            v_ind = np.zeros([1, r.shape[-1]], dtype = complex)
            n_ind = np.zeros([1, r.shape[-1]], dtype = complex)
        elif sum == '1D':            
            self.printtxt('1D sum not implemented')
            return 
        else:
            v_ind = np.zeros([1, r.shape[1], r.shape[2], r.shape[3]], dtype = complex)
            n_ind = np.zeros([1, r.shape[1], r.shape[2], r.shape[3]], dtype = complex)

        dm_wGG = self.dm_wGG
        eps_GG_plus = dm_wGG[0]
        eig_plus, vec_plus = np.linalg.eig(eps_GG_plus)  # find eigenvalues and eigenvectors
        vec_plus_dual = np.linalg.inv(vec_plus)
        
        # loop over frequencies, where the eigenvalues for the 2D matrix in G,G' are found.
        for i in np.array(range(self.Nw-1))+1:
            eps_GG = eps_GG_plus
            eig, vec = eig_plus,vec_plus
            vec_dual = vec_plus_dual
            eps_GG_plus = dm_wGG[i] # epsilon_GG'(omega + d-omega)
            eig_plus, vec_plus = np.linalg.eig(eps_GG_plus)
            vec_plus_dual = np.linalg.inv(vec_plus)
            eig_dummy[0,:] = eig
            eig_all = np.append(eig_all, eig_dummy, axis=0) # append all eigenvalues to array         
            # loop to check find the eigenvalues that crosses zero from negative to positive values:
            for k in range(self.npw):
                for m in range(self.npw):
                    if eig[k]< 0 and 0 < eig_plus[m]:
                        # check it's the same mode - Overlap between eigenvectors should be large:
                        if abs(np.inner(vec[:,k], np.conjugate(vec_plus[:,m]))) > 0.95:                             
                            self.printtxt('crossing found at w = %1.1f eV'%self.w_w[i-1])
                            eig_left = np.append(eig_left, eig[k])   
                            eig_right = np.append(eig_right, eig_plus[m])

                            vec_dummy[0, :] = vec[:,k]
                            vec_modes = np.append(vec_modes, vec_dummy, axis = 0)
                            vec_dummy[0, :] = vec_dual[k, :].T
                            vec_modes_dual = np.append(vec_modes_dual, vec_dummy, axis = 0)
                            
                            w1 = self.w_w[i-1]
                            w2 = self.w_w[i]
                            a = np.real((eig_plus[m]-eig[k]) / (w2-w1))
                            w0 = np.real(-eig[k]) / a + w1
                            eig0 = a*(w0-w1)+eig[k]

                            w_0 = np.append(w_0,w0)
                            w_left = np.append(w_left, w1)
                            eig_0 = np.append(eig_0,eig0)                           

                            n_dummy = np.zeros([1, r.shape[1], r.shape[2],
                                                r.shape[3]], dtype = complex)
                            v_dummy = np.zeros([1, r.shape[1], r.shape[2],
                                                r.shape[3]], dtype = complex)    

                            vec_n = np.zeros([self.npw])
                            for iG in range(self.npw):  # Fourier transform
                                qG = np.dot((q + self.Gvec_Gc[iG]), self.bcell_cv)
                                coef_G = np.dot(qG, qG) / (4 * pi)
                                qGr_R = np.inner(qG, r.T).T
                                v_dummy += vec[iG, k] * np.exp(1j * qGr_R)
                                n_dummy += vec[iG, k] * np.exp(1j * qGr_R) * coef_G

                            if checkphase:
                                integral = np.zeros([81])
                                integral_n = np.zeros([81])
                                phases = np.linspace(0,2,81)
                                phases_n = np.linspace(0,2,81)
                                for ip in range(81):
                                    v_int = v_dummy * np.exp(1j * pi * phases[ip])
                                    integral[ip] = abs(np.imag(v_int)).sum()
                                    n_int = n_dummy * np.exp(1j * pi * phases_n[ip])
                                    integral_n[ip] = abs(np.imag(n_int)).sum()
                                    
                                phase = phases[np.argsort(integral)][0]
                                phase_n = phases_n[np.argsort(integral_n)][0]
                                v_dummy = v_dummy * np.exp(1j * pi * phase)
                                n_dummy = n_dummy * np.exp(1j * pi * phase_n)
#
                            if sum == '2D':
                                i_xyz = 3 
                                v_dummy_z = np.zeros([1,v_dummy.shape[i_xyz]],
                                                     dtype = complex)
                                n_dummy_z = np.zeros([1,v_dummy.shape[i_xyz]],
                                                     dtype = complex)
                                v_dummy_z[0,:] = np.sum(np.sum(v_dummy, axis = 1),
                                                        axis = 1)[0,:]
                                n_dummy_z[0,:] = np.sum(np.sum(n_dummy, axis = 1),
                                                        axis = 1)[0,:]

                                v_ind = np.append(v_ind, v_dummy_z, axis=0)
                                n_ind = np.append(n_ind, n_dummy_z, axis=0)
                                
                            #elif sum == '1D':
                            #    self.printtxt('ALDA not implemented')
                                
                            elif sum == '1D':
                                self.printtxt('1D sum not implemented')
                            else :
                                v_ind = np.append(v_ind, v_dummy, axis=0)
                                n_ind = np.append(n_ind, n_dummy, axis=0)                            
        """                        
        returns: grid points, frequency grid, all eigenvalues, mode energies,
                 mode eigenvalues, eigenvalues of left and right-side points,
                 mode eigenvectors, mode dual eigenvectors,
                 induced potential in real space, induced density in real space
        """
        return rrr, self.w_w, eig_all[1:], w_0, eig_0, w_left, eig_left, \
               eig_right, vec_modes[1:], vec_modes_dual[1:], v_ind[1:], n_ind[1:]
    

    def project_chi_to_LCAO_pair_orbital(self, orb_MG):

        nLCAO = orb_MG.shape[0]
        N = np.zeros((self.Nw, nLCAO, nLCAO), dtype=complex)

        kcoulinv_GG = np.zeros((self.npw, self.npw))
        for iG in range(self.npw):
            qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
            kcoulinv_GG[iG, iG] = np.dot(qG, qG)

        kcoulinv_GG /= 4.*pi

        dm_wGG = self.get_RPA_dielectric_matrix()

        for mu in range(nLCAO):
            for nu in range(nLCAO):
                pairorb_R = orb_MG[mu] * orb_MG[nu]
                if not (pairorb_R * pairorb_R.conj() < 1e-10).all():
                    tmp_G = np.fft.fftn(pairorb_R) * self.vol / self.nG0

                    pairorb_G = np.zeros(self.npw, dtype=complex)
                    for iG in range(self.npw):
                        index = self.Gindex[iG]
                        pairorb_G[iG] = tmp_G[index[0], index[1], index[2]]

                    for iw in range(self.Nw):
                        chi_GG = (dm_wGG[iw] - np.eye(self.npw)) * kcoulinv_GG
                        N[iw, mu, nu] = (np.outer(pairorb_G.conj(), pairorb_G) * chi_GG).sum()
#                        N[iw, mu, nu] = np.inner(pairorb_G.conj(),np.inner(pairorb_G, chi_GG))

        return N


    def write(self, filename, all=False):
        """Dump essential data"""

        data = {'nbands': self.nbands,
                'acell': self.acell_cv, #* Bohr,
                'bcell': self.bcell_cv, #/ Bohr,
                'h_cv' : self.h_cv,   #* Bohr,
                'nG'   : self.nG,
                'nG0'  : self.nG0,
                'vol'  : self.vol,   #* Bohr**3,
                'BZvol': self.BZvol, #/ Bohr**3,
                'nkpt' : self.nkpt,
                'ecut' : self.ecut,  #* Hartree,
                'npw'  : self.npw,
                'eta'  : self.eta,   #* Hartree,
                'ftol' : self.ftol,  #* self.nkpt,
                'Nw'   : self.Nw,
                'NwS'  : self.NwS,
                'dw'   : self.dw,    # * Hartree,
                'q_red': self.q_c,
                'q_car': self.qq_v,    # / Bohr,
                'qmod' : np.dot(self.qq_v, self.qq_v), # / Bohr
                #'vcut' : self.vcut,
                'nvalence'     : self.nvalence,                
                'hilbert_trans' : self.hilbert_trans,
                'optical_limit' : self.optical_limit,
                'e_skn'         : self.e_skn,          # * Hartree,
                'f_skn'         : self.f_skn,          # * self.nkpt,
                'bzk_kc'       : self.bzk_kc,
                'ibzk_kc'      : self.ibzk_kc,
                'kq_k'         : self.kq_k,
                'Gvec_Gc'      : self.Gvec_Gc,
                'dfNLFRPA_w'   : self.df1_w,
                'dfLFCRPA_w'   : self.df2_w,
                'dfNLFALDA_w'  : self.df3_w,
                'dfLFCALDA_w'  : self.df4_w,
                'df_flag'      : True}

        if all == True:
            from gpaw.response.parallel import par_write
            par_write('chi0' + filename,'chi0_wGG',self.wcomm,self.chi0_wGG)
        
        if rank == 0:
            pickle.dump(data, open(filename, 'w'), -1)

        self.comm.barrier()


    def read(self, filename):
        """Read data from pickle file"""

        data = pickle.load(open(filename))
        
        self.nbands = data['nbands']
        self.acell_cv = data['acell']
        self.bcell_cv = data['bcell']
        self.h_cv   = data['h_cv']
        self.nG    = data['nG']
        self.nG0   = data['nG0']
        self.vol   = data['vol']
        self.BZvol = data['BZvol']
        self.nkpt  = data['nkpt']
        self.ecut  = data['ecut']
        self.npw   = data['npw']
        self.eta   = data['eta']
        self.ftol  = data['ftol']
        self.Nw    = data['Nw']
        self.NwS   = data['NwS']
        self.dw    = data['dw']
        self.q_c   = data['q_red']
        self.qq_v  = data['q_car']
        self.qmod  = data['qmod']
        #self.vcut  = data['vcut']
        
        self.hilbert_trans = data['hilbert_trans']
        self.optical_limit = data['optical_limit']
        self.e_skn  = data['e_skn']
        self.f_skn  = data['f_skn']
        self.nvalence= data['nvalence']
        self.bzk_kc  = data['bzk_kc']
        self.ibzk_kc = data['ibzk_kc']
        self.kq_k    = data['kq_k']
        self.Gvec_Gc  = data['Gvec_Gc']
        self.df1_w   = data['dfNLFRPA_w']
        self.df2_w   = data['dfLFCRPA_w']
        self.df3_w   = data['dfNLFALDA_w']
        self.df4_w   = data['dfLFCALDA_w']
        self.df_flag = data['df_flag']
        
        self.printtxt('Read succesfully !')
