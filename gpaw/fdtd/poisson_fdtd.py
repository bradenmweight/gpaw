# Electrodynamics module, by Arto Sakko (Aalto University)

from ase.parallel import parprint
from ase.units import Hartree, Bohr, _eps0, _c, _aut
from gpaw import PoissonConvergenceError
from gpaw.grid_descriptor import GridDescriptor
from gpaw.io import open as gpaw_io_open
from gpaw.mpi import world, serial_comm
from gpaw.tddft.units import attosec_to_autime, autime_to_attosec
from gpaw.transformers import Transformer
from gpaw.utilities.blas import axpy
from gpaw.utilities.gpts import get_number_of_grid_points
from math import pi
from gpaw.utilities.gauss import Gaussian
from poisson_corr import PoissonSolver
from polarizable_material import *
from potential_couplers import *
from string import split
import _gpaw
import gpaw.mpi as mpi
import numpy as np
import sys


# In atomic units, 1/(4*pi*e_0) = 1
_maxL = 1  # 1 for monopole, 4 for dipole, 9 for quadrupole

# This helps in telling the classical quantitites from the quantum ones
class PoissonOrganizer:
    def __init__(self, poisson_solver=None):
        self.poisson_solver = poisson_solver
        self.gd = None
        self.density = None
        self.cell = None
        self.spacing_def = None
        self.spacing = None

# Contains one PoissonSolver for the classical and one for the quantum subsystem
class FDTDPoissonSolver:
    def __init__(self, nn=3,
                        relax='J',
                        eps=2e-10,
                        classical_material=None,
                        cell=None,
                        qm_spacing=0.30,
                        cl_spacing=1.20,
                        tag='fdtd.poisson',
                        remove_moments=(_maxL, 1),
                        potential_coupler='Refiner',
                        coupling_level='both',
                        communicator=serial_comm,
                        debug_plots=0):

        assert(coupling_level in ['none', 'both', 'cl2qm', 'qm2cl'])
        self.potential_coupling_scheme = potential_coupler
        self.coupling_level = coupling_level
        
        if classical_material == None:
            self.classical_material = PolarizableMaterial()
        else:
            self.classical_material = classical_material
        
        self.description = 'FDTD+TDDFT'
        self.set_calculation_mode('solve')
        
        self.remove_moment_qm = remove_moments[0]
        self.remove_moment_cl = remove_moments[1]
        self.tag = tag
        self.time = 0.0
        self.time_step = 0.0
        self.rank = mpi.rank
        self.dm_file = None
        self.kick = None
        self.debug_plots = debug_plots
        self.maxiter = 2000
        
        # Only handle the quantities via self.qm or self.cl 
        self.cl = PoissonOrganizer()        
        self.cl.spacing_def = cl_spacing * np.ones(3) / Bohr
        self.cl.extrapolated_qm_phi = None
        self.cl.dcomm = communicator
        self.cl.dparsize = None
        self.qm = PoissonOrganizer(PoissonSolver)  # Default solver
        self.qm.spacing_def = qm_spacing * np.ones(3) / Bohr
        self.qm.cell = np.array(cell) / Bohr
        
        # Create grid descriptor for the classical part
        _cell = np.array(cell) / Bohr
        self.cl.spacing = self.cl.spacing_def
        if np.size(_cell) == 3:
            self.cl.cell = np.diag(_cell)
        else:
            self.cl.cell = _cell

        N_c = get_number_of_grid_points(self.cl.cell, self.cl.spacing)
        self.cl.spacing = np.diag(self.cl.cell) / N_c
        self.cl.gd = GridDescriptor(N_c,
                                    self.cl.cell,
                                    False,
                                    self.cl.dcomm,
                                    self.cl.dparsize)
        self.cl.gd_global = GridDescriptor(N_c,
                                           self.cl.cell,
                                           False,
                                           serial_comm,
                                           None)
        self.cl.extrapolated_qm_phi = self.cl.gd.empty()
        
        parprint('FDTDPoissonSolver: domain parallelization with %i processes.' % self.cl.gd.comm.size)

    def estimate_memory(self, mem):
        #self.cl.poisson_solver.estimate_memory(mem)
        self.qm.poisson_solver.estimate_memory(mem)

    # Return the TDDFT stencil by default 
    def get_stencil(self, mode='qm'):
        if mode=='qm':
            return self.qm.poisson_solver.get_stencil()
        else:
            return self.cl.poisson_solver.get_stencil()

    # Initialize both PoissonSolvers
    def initialize(self, load_Gauss=False):
        self.qm.poisson_solver.initialize(load_Gauss)
        self.cl.poisson_solver.initialize(load_Gauss)     

    def set_grid_descriptor(self, qmgd):

        self.qm.gd = qmgd
        
        # Create quantum Poisson solver
        self.qm.poisson_solver = PoissonSolver()
        self.qm.poisson_solver.set_grid_descriptor(self.qm.gd)
        self.qm.poisson_solver.initialize()
        self.qm.phi = self.qm.gd.zeros()
        self.qm.rho = self.qm.gd.zeros()

        # Set quantum grid descriptor
        self.qm.poisson_solver.set_grid_descriptor(qmgd)

        # Create classical PoissonSolver
        self.cl.poisson_solver = PoissonSolver()
        self.cl.poisson_solver.set_grid_descriptor(self.cl.gd)
        self.cl.poisson_solver.initialize()
            
        # Initialize classical material, its Poisson solver was generated already
        self.cl.poisson_solver.set_grid_descriptor(self.cl.gd)
        self.classical_material.initialize(self.cl.gd)
        self.cl.extrapolated_qm_phi = self.cl.gd.zeros()
        self.cl.phi = self.cl.gd.zeros()
        self.cl.extrapolated_qm_phi = self.cl.gd.empty()

        # Initialize potential coupler
        if self.potential_coupling_scheme == 'Multipoles':
            parprint('Classical-quantum coupling by multipole expansion with maxL: %i and coupling level: %s' % (self.remove_moment_qm, self.coupling_level))
            self.potential_coupler = MultipolesPotentialCoupler(qm = self.qm,
                                                                cl = self.cl,
                                                                index_offset_1 = self.shift_indices_1,
                                                                index_offset_2 = self.shift_indices_2,
                                                                extended_index_offset_1 = self.extended_shift_indices_1,
                                                                extended_index_offset_2 = self.extended_shift_indices_2,
                                                                extended_delta_index = self.extended_deltaIndex,
                                                                num_refinements = self.num_refinements,
                                                                remove_moment_qm = self.remove_moment_qm,
                                                                remove_moment_cl = self.remove_moment_cl,
                                                                coupling_level = self.coupling_level,
                                                                rank = self.rank)
        else:
            parprint('Classical-quantum coupling by coarsening/refining')
            self.potential_coupler = RefinerPotentialCoupler(qm = self.qm,
                                                             cl = self.cl,
                                                             index_offset_1 = self.shift_indices_1,
                                                             index_offset_2 = self.shift_indices_2,
                                                             extended_index_offset_1 = self.extended_shift_indices_1,
                                                             extended_index_offset_2 = self.extended_shift_indices_2,
                                                             extended_delta_index = self.extended_deltaIndex,
                                                             num_refinements = self.num_refinements,
                                                             remove_moment_qm = self.remove_moment_qm,
                                                             remove_moment_cl = self.remove_moment_cl,
                                                             coupling_level = self.coupling_level,
                                                             rank = self.rank)
            
        self.phi_tot_clgd = self.cl.gd.empty()
        self.phi_tot_qmgd = self.qm.gd.empty()

    def cut_cell(self, atoms_in, vacuum=5.0, corners=None):
        qmh = self.qm.spacing_def
        if corners != None:
            v1 = np.array(corners[0]).ravel() / Bohr
            v2 = np.array(corners[1]).ravel() / Bohr
        else: # Use vacuum
            pos_old = atoms_in.get_positions()[0];
            dmy_atoms = atoms_in.copy()
            dmy_atoms.center(vacuum=vacuum)
            pos_new = dmy_atoms.get_positions()[0];
            v1 = (pos_old - pos_new)/Bohr
            v2 = v1 + np.diag(dmy_atoms.get_cell())/Bohr

        # Sanity check: quantum box must be inside the classical one
        assert(all([v1[w] <= v2[w] and
                    v1[w] >= 0 and
                    v2[w] <= np.diag(self.cl.cell)[w] for w in range(3)]))
        
        # Create new Atoms object
        atoms_out = atoms_in.copy()

        # Quantum grid is probably not yet created
        if not self.qm.gd:
            self.qm.cell = np.zeros((3, 3))
            for w in range(3):
                self.qm.cell[w, w] = v2[w] - v1[w]
        
            N_c = get_number_of_grid_points(self.qm.cell, qmh)
            self.qm.spacing = np.diag(self.qm.cell) / N_c
        else:
            self.qm.cell = self.qm.gd.cell_cv
            N_c = self.qm.gd.N_c
            self.qm.spacing = self.qm.gd.get_grid_spacings()
        
        # Ratios of the user-given spacings
        hratios = self.cl.spacing_def / qmh
        self.num_refinements = 1 + int(round(np.log(hratios[0]) / np.log(2.0)))
        assert([int(round(np.log(hratios[w]) / np.log(2.0))) == self.num_refinements for w in range(3)])

        # Classical corner indices must be divisable with numb
        if any(self.cl.spacing / self.qm.spacing >= 3):
            numb = 1
        elif any(self.cl.spacing / self.qm.spacing >= 2):
            numb = 2
        else:
            numb = 4
        
        # The index mismatch of the two simulation cells
        self.num_indices = numb * np.ceil((np.array(v2) -
                                           np.array(v1)) /
                                          self.cl.spacing / numb)
        
        self.num_indices_1 = numb * np.floor(np.array(v1) / self.cl.spacing / numb)
        self.num_indices_2 = numb * np.ceil(np.array(v2) / self.cl.spacing / numb)
        self.num_indices = self.num_indices_2 - self.num_indices_1
        
        # Center, left, and right points of the suggested quantum grid
        cp = 0.5 * (np.array(v1) + np.array(v2))
        lp = cp - 0.5 * self.num_indices * self.cl.spacing 
        rp = cp + 0.5 * self.num_indices * self.cl.spacing
                
        # Indices in the classical grid restricting the quantum grid
        self.shift_indices_1 = np.round(lp / self.cl.spacing)
        self.shift_indices_2 = self.shift_indices_1 + self.num_indices

        # Sanity checks
        assert(all([self.shift_indices_1[w] >= 0 and
                    self.shift_indices_2[w] <= self.cl.gd.N_c[w] for w in range(3)])), \
                    "Could not find appropriate quantum grid. Move it further away from the boundary."
        
        # Corner coordinates
        self.qm.corner1 = self.shift_indices_1 * self.cl.spacing
        self.qm.corner2 = self.shift_indices_2 * self.cl.spacing
        
        # New quantum grid
        for w in range(3):
            self.qm.cell[w, w] = (self.shift_indices_2[w] -
                                  self.shift_indices_1[w]) * \
                                  self.cl.spacing[w]
        self.qm.spacing = self.cl.spacing / hratios
        N_c = get_number_of_grid_points(self.qm.cell, self.qm.spacing)
        
        atoms_out.set_cell(np.diag(self.qm.cell) * Bohr)
        atoms_out.positions = atoms_in.get_positions() - self.qm.corner1 * Bohr
        
        parprint("Quantum box readjustment:")
        parprint("  Given cell/atomic coordinates:")
        parprint("             [%10.5f %10.5f %10.5f]" % tuple(np.diag(atoms_in.get_cell())))
        for s, c in zip(atoms_in.get_chemical_symbols(), atoms_in.get_positions()):
            parprint("           %s %10.5f %10.5f %10.5f" % (s, c[0], c[1], c[2]))
        parprint("  Readjusted cell/atomic coordinates:")
        parprint("             [%10.5f %10.5f %10.5f]" % tuple(np.diag(atoms_out.get_cell())))
        for s, c in zip(atoms_out.get_chemical_symbols(), atoms_out.get_positions()):
            parprint("           %s %10.5f %10.5f %10.5f" % (s, c[0], c[1], c[2]))
        
        parprint("  Given corner points:       (%10.5f %10.5f %10.5f) - (%10.5f %10.5f %10.5f)" %
                 (tuple(np.concatenate((v1, v2)) * Bohr)))
        parprint("  Readjusted corner points:  (%10.5f %10.5f %10.5f) - (%10.5f %10.5f %10.5f)" %
                 (tuple(np.concatenate((self.qm.corner1,
                                        self.qm.corner2)) * Bohr)))
        parprint("  Indices in classical grid: (%10i %10i %10i) - (%10i %10i %10i)" %
                 (tuple(np.concatenate((self.shift_indices_1,
                                        self.shift_indices_2)))))
        parprint("  Grid points in classical grid: (%10i %10i %10i)" % (tuple(self.cl.gd.N_c)))
        parprint("  Grid points in quantum grid:   (%10i %10i %10i)" % (tuple(N_c)))
        
        parprint("  Spacings in quantum grid:    (%10.5f %10.5f %10.5f)" %
                 (tuple(np.diag(self.qm.cell) * Bohr / N_c)))
        parprint("  Spacings in classical grid:  (%10.5f %10.5f %10.5f)" %
                 (tuple(np.diag(self.cl.cell) * Bohr / \
                        get_number_of_grid_points(self.cl.cell, self.cl.spacing))))
        parprint("  Ratios of cl/qm spacings:    (%10i %10i %10i)" % (tuple(hratios)))
        parprint("                             = (%10.2f %10.2f %10.2f)" %
                 (tuple((np.diag(self.cl.cell) * Bohr / \
                         get_number_of_grid_points(self.cl.cell,
                                                   self.cl.spacing)) / \
                        (np.diag(self.qm.cell) * Bohr / N_c))))
        parprint("  Needed number of refinements: %10i" % self.num_refinements)
        
        #   First, create the quantum grid equivalent griddescriptor object self.cl.subgd.
        #   Then coarsen it until its h_cv equals that of self.cl.gd.
        #   Finally, map the points between clgd and coarsened subgrid.
        subcell_cv = np.diag(self.qm.corner2 - self.qm.corner1)
        N_c = get_number_of_grid_points(subcell_cv, self.cl.spacing)
        N_c = self.shift_indices_2 - self.shift_indices_1
        self.cl.subgds = []
        self.cl.subgds.append(GridDescriptor(N_c, subcell_cv, False, serial_comm, self.cl.dparsize))

        parprint("  N_c/spacing of the subgrid:           %3i %3i %3i / %.4f %.4f %.4f" % 
                  (self.cl.subgds[0].N_c[0],
                   self.cl.subgds[0].N_c[1],
                   self.cl.subgds[0].N_c[2],
                   self.cl.subgds[0].h_cv[0][0] * Bohr,
                   self.cl.subgds[0].h_cv[1][1] * Bohr,
                   self.cl.subgds[0].h_cv[2][2] * Bohr))
        parprint("  shape from the subgrid:           %3i %3i %3i" % (tuple(self.cl.subgds[0].empty().shape)))

        self.cl.coarseners = []
        self.cl.refiners = []
        for n in range(self.num_refinements):
            self.cl.subgds.append(self.cl.subgds[n].refine())
            self.cl.refiners.append(Transformer(self.cl.subgds[n], self.cl.subgds[n + 1]))
            
            parprint("  refiners[%i] can perform the transformation (%3i %3i %3i) -> (%3i %3i %3i)" % (\
                     n,
                     self.cl.subgds[n].empty().shape[0],
                     self.cl.subgds[n].empty().shape[1],
                     self.cl.subgds[n].empty().shape[2],
                     self.cl.subgds[n + 1].empty().shape[0],
                     self.cl.subgds[n + 1].empty().shape[1],
                     self.cl.subgds[n + 1].empty().shape[2]))
            self.cl.coarseners.append(Transformer(self.cl.subgds[n + 1], self.cl.subgds[n]))
        self.cl.coarseners[:] = self.cl.coarseners[::-1]
        
        # Now extend the grid in order to handle the zero boundary conditions that the refiner assumes
        # The default interpolation order
        self.extend_nn = Transformer(GridDescriptor([8, 8, 8], [1, 1, 1], False, serial_comm, None),
                                     GridDescriptor([8, 8, 8], [1, 1, 1], False, serial_comm, None).coarsen()).nn
        
        self.extended_num_indices = self.num_indices + [2, 2, 2]
        
        # Center, left, and right points of the suggested quantum grid
        extended_cp = 0.5 * (np.array(v1) + np.array(v2))
        extended_lp = extended_cp - 0.5 * (self.extended_num_indices) * self.cl.spacing 
        extended_rp = extended_cp + 0.5 * (self.extended_num_indices) * self.cl.spacing
        
        # Indices in the classical grid restricting the quantum grid
        self.extended_shift_indices_1 = np.floor(extended_lp / self.cl.spacing)
        self.extended_shift_indices_2 = self.extended_shift_indices_1 + self.extended_num_indices

        # Sanity checks
        assert(all([self.extended_shift_indices_1[w] >= 0 and
                    self.extended_shift_indices_2[w] <= self.cl.gd.N_c[w] for w in range(3)])), \
                    "Could not find appropriate quantum grid. Move it further away from the boundary."
        
        # Corner coordinates
        self.qm.extended_corner1 = self.extended_shift_indices_1 * self.cl.spacing
        self.qm.extended_corner2 = self.extended_shift_indices_2 * self.cl.spacing
        N_c = self.extended_shift_indices_2 - self.extended_shift_indices_1
               
        self.cl.extended_subgds = []
        self.cl.extended_refiners = []
        extended_subcell_cv = np.diag(self.qm.extended_corner2 - self.qm.extended_corner1)

        self.cl.extended_subgds.append(GridDescriptor(N_c,
                                                      extended_subcell_cv,
                                                      False,
                                                      serial_comm,
                                                      None))
        
        for n in range(self.num_refinements):
            self.cl.extended_subgds.append(self.cl.extended_subgds[n].refine())
            self.cl.extended_refiners.append(Transformer(self.cl.extended_subgds[n], self.cl.extended_subgds[n + 1]))
            parprint("  extended_refiners[%i] can perform the transformation (%3i %3i %3i) -> (%3i %3i %3i)" %
                    (n,
                     self.cl.extended_subgds[n].empty().shape[0],
                     self.cl.extended_subgds[n].empty().shape[1],
                     self.cl.extended_subgds[n].empty().shape[2],
                     self.cl.extended_subgds[n + 1].empty().shape[0],
                     self.cl.extended_subgds[n + 1].empty().shape[1],
                     self.cl.extended_subgds[n + 1].empty().shape[2]))
        
        parprint("  N_c/spacing of the refined subgrid:   %3i %3i %3i / %.4f %.4f %.4f" % 
                  (self.cl.subgds[-1].N_c[0],
                   self.cl.subgds[-1].N_c[1],
                   self.cl.subgds[-1].N_c[2],
                   self.cl.subgds[-1].h_cv[0][0] * Bohr,
                   self.cl.subgds[-1].h_cv[1][1] * Bohr,
                   self.cl.subgds[-1].h_cv[2][2] * Bohr))
        parprint("  shape from the refined subgrid:       %3i %3i %3i" % 
                 (tuple(self.cl.subgds[-1].empty().shape)))
        
        self.extended_deltaIndex = 2 ** (self.num_refinements) * self.extend_nn
        parprint("self.extended_deltaIndex = %i" % self.extended_deltaIndex)
        
        qgpts = self.cl.subgds[-2].N_c
        
        # Assure that one returns to the original shape
        dmygd = self.cl.subgds[-1].coarsen()
        for n in range(self.num_refinements - 1):
            dmygd = dmygd.coarsen()
        
        parprint("  N_c/spacing of the coarsened subgrid: %3i %3i %3i / %.4f %.4f %.4f" % 
                  (dmygd.N_c[0], dmygd.N_c[1], dmygd.N_c[2],
                   dmygd.h_cv[0][0] * Bohr, dmygd.h_cv[1][1] * Bohr, dmygd.h_cv[2][2] * Bohr))
       
        return atoms_out, self.qm.spacing[0] * Bohr, qgpts

   
    # Where the induced dipole moment is written
    def set_dipole_moment_fname(self, fname):
        self.fname = fname

    # Set the time step
    def set_time_step(self, time_step):
        self.time_step = time_step

    # This must be called before propagation begins
    def initialize_propagation(self, kick,
                                     time=0.0):
        self.time = time
        self.kick = kick
        
        # dipole moment file
        if self.rank == 0:
            if self.dm_file is not None and not self.dm_file.closed:
                raise RuntimeError('Dipole moment file is already open')
            if self.time == 0.0:
                mode = 'w'
            else:
                mode = 'a'
            self.dm_file = file(self.fname, mode)
            if self.dm_file.tell() == 0:
                header = '# Kick = [%22.12le, %22.12le, %22.12le]\n' % \
                            (self.kick[0], self.kick[1], self.kick[2])
                header += '# %15s %15s %22s %22s %22s\n' % \
                            ('time', 'norm', 'dmx', 'dmy', 'dmz')
                self.dm_file.write(header)
                self.dm_file.flush()

    def finalize_propagation(self):
        if self.rank == 0:
            self.dm_file.close()
            self.dm_file = None
    
    def set_calculation_mode(self, calculation_mode):
        # Three calculation modes are available:
        #  1) solve:     just solve the Poisson equation with
        #                given quantum+classical rho
        #  2) iterate:   iterate classical density so that the Poisson
        #                equation for quantum+classical rho is satisfied
        #  3) propagate: propagate classical density in time, and solve
        #                the new Poisson equation
        assert(calculation_mode == 'solve' or
               calculation_mode == 'iterate' or
               calculation_mode == 'propagate')
        self.calculation_mode = calculation_mode

    # The density object must be attached, so that the electric field
    # from all-electron density can be calculated    
    def set_density(self, density):
        self.density = density
        
    # Returns the classical density and the grid descriptor
    def get_density(self, global_array=False):
        if global_array:
            return self.cl.gd.collect(self.classical_material.charge_density) * \
                   self.classical_material.sign, \
                   self.cl.gd
        else:
            return self.classical_material.charge_density * \
                   self.classical_material.sign, \
                   self.cl.gd
        
    # Returns the quantum + classical density in the large classical box,
    # so that the classical charge is coarsened into it and the quantum
    # charge is refined there
    def get_combined_data(self, qmdata=None, cldata=None, spacing=None):
        
        if qmdata == None:
            qmdata = self.density.rhot_g
        
        if cldata == None:
            cldata = self.classical_material.charge_density
        
        if spacing == None:
            spacing = self.cl.gd.h_cv[0, 0]
        
        spacing_au = spacing / Bohr  # from Angstroms to a.u.
        
        # Collect data from different processes
        cln = self.cl.gd.collect(cldata) * self.classical_material.sign
        qmn = self.qm.gd.collect(qmdata)

        clgd = GridDescriptor(self.cl.gd.N_c,
                              self.cl.cell,
                              False,
                              serial_comm,
                              None)

        if world.rank == 0:
            # refine classical part
            while clgd.h_cv[0, 0] > spacing_au * 1.50:  # 45:
                cln = Transformer(clgd, clgd.refine()).apply(cln)
                clgd = clgd.refine()
                
            # refine quantum part
            qmgd = GridDescriptor(self.qm.gd.N_c,
                                  self.qm.cell,
                                  False,
                                  serial_comm,
                                  None)                           
            while qmgd.h_cv[0, 0] < clgd.h_cv[0, 0] * 0.95:
                qmn = Transformer(qmgd, qmgd.coarsen()).apply(qmn)
                qmgd = qmgd.coarsen()
            
            assert np.all(qmgd.h_cv == clgd.h_cv), " Spacings %.8f (qm) and %.8f (cl) Angstroms" % (qmgd.h_cv[0][0] * Bohr, clgd.h_cv[0][0] * Bohr)
            
            # now find the corners
            r_gv_cl = clgd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
            cind = self.qm.corner1 / np.diag(clgd.h_cv) - 1
            
            n = qmn.shape

            # print 'Corner points:     ', self.qm.corner1*Bohr,      ' - ', self.qm.corner2*Bohr
            # print 'Calculated points: ', r_gv_cl[tuple(cind)]*Bohr, ' - ', r_gv_cl[tuple(cind+n+1)]*Bohr
                        
            cln[cind[0] + 1:cind[0] + n[0] + 1,
                cind[1] + 1:cind[1] + n[1] + 1,
                cind[2] + 1:cind[2] + n[2] + 1] += qmn
        
        world.barrier()
        return cln, clgd
            
    
    # Solve quantum and classical potentials, and add them up
    def solve_solve(self, **kwargs):
        self.phi_tot_qmgd, self.phi_tot_clgd, niter = self.potential_coupler.getPotential(local_rho_qm_qmgd = self.qm.rho, local_rho_cl_clgd = self.classical_material.sign * self.classical_material.charge_density, **kwargs)
        self.qm.phi[:] = self.phi_tot_qmgd[:]
        self.cl.phi[:] = self.phi_tot_clgd[:]
        return niter

 
    # Iterate classical and quantum potentials until convergence
    def solve_iterate(self, **kwargs):
        # Initial value (unefficient?) 
        self.solve_solve(**kwargs)
        old_rho_qm = self.qm.rho.copy()
        old_rho_cl = self.classical_material.charge_density.copy()
        
        niter_cl = 0
            
        while True:
            # field from the potential
            self.classical_material.solve_electric_field(self.cl.phi)  # E = -Div[Vh]

            # Polarizations P0_j and Ptot
            self.classical_material.solve_polarizations()  # P = (eps - eps0)E

            # Classical charge density
            self.classical_material.solve_rho()  # n = -Grad[P]
                
            # Update electrostatic potential         # nabla^2 Vh = -4*pi*n
            niter = self.solve_solve(**kwargs) 

            # Mix potential
            try:
                self.mix_phi
            except:
                self.mix_phi = SimpleMixer(0.10, self.qm.phi)

            self.qm.phi = self.mix_phi.mix(self.qm.phi)
                
            # Check convergence
            niter_cl += 1
            
            dRho = self.qm.gd.integrate(abs(self.qm.rho - old_rho_qm)) + \
                    self.cl.gd.integrate(abs(self.classical_material.charge_density - old_rho_cl))
            
            if(abs(dRho) < 1e-3):
                break
            old_rho_qm = rho.copy()
            old_rho_cl = (self.classical_material.sign * self.classical_material.charge_density).copy()

        return (niter, niter_cl)
        
            
    def solve_propagate(self, **kwargs):        
        if self.debug_plots!=0 and np.floor(self.time / self.time_step) % self.debug_plots == 0:
            from visualization import visualize_density
            visualize_density(self, plotInduced=False)

        # 1) P(t) from P(t-dt) and J(t-dt/2)
        self.classical_material.propagate_polarizations(self.time_step)
                
        # 2) n(t) from P(t)
        self.classical_material.solve_rho()
        
        # 3a) V(t) from n(t)
        niter = self.solve_solve(**kwargs)
        
        # 4a) E(r) from V(t):      E = -Div[Vh]
        self.classical_material.solve_electric_field(self.cl.phi)
                
        # 4b) Apply the kick by changing the electric field
        if self.time == 0:
            self.classical_material.kick_electric_field(self.time_step, self.kick)
                    
        # 5) J(t+dt/2) from J(t-dt/2) and P(t)
        self.classical_material.propagate_currents(self.time_step)

        # Write updated dipole moment into file
        self.update_dipole_moment_file(self.qm.rho)
                
        # Update timer
        self.time = self.time + self.time_step
                
        # Do not propagate before the next time step
        self.set_calculation_mode('solve')

        return niter
                                

    def solve(self, phi,
                    rho,
                    charge=None,
                    eps=None,
                    maxcharge=1e-6,
                    zero_initial_phi=False,
                    calculation_mode=None):

        if self.density == None:
            print 'FDTDPoissonSolver requires a density object.' \
                  ' Use set_density routine to initialize it.'
            raise

        # Update local variables (which may have changed in SCF cycle or propagator) 
        self.qm.phi = phi
        self.qm.rho = rho

        if(self.calculation_mode == 'solve'):  # do not modify the polarizable material
            return self.solve_solve(charge=None,
                                   eps=None,
                                   maxcharge=maxcharge,
                                   zero_initial_phi=False)

        elif(self.calculation_mode == 'iterate'):  # find self-consistent density
            return self.solve_iterate(charge=None,
                                      eps=None,
                                      maxcharge=maxcharge,
                                      zero_initial_phi=False)
        
        elif(self.calculation_mode == 'propagate'):  # propagate one time step
            return self.solve_propagate(charge=None,
                                        eps=None,
                                        maxcharge=maxcharge,
                                        zero_initial_phi=False)

        
    def update_dipole_moment_file(self, rho):
        # Classical contribution. Note the different origin.
        r_gv = self.cl.gd.get_grid_point_coordinates().transpose((1, 2, 3, 0)) - self.qm.corner1
        dmcl = -1.0 * self.classical_material.sign * np.array([self.cl.gd.integrate(np.multiply(r_gv[:, :, :, w] + self.qm.corner1[w], self.classical_material.charge_density)) for w in range(3)])

        # Quantum contribution
        dm = self.density.finegd.calculate_dipole_moment(self.density.rhot_g) + dmcl
        norm = self.qm.gd.integrate(rho) + self.classical_material.sign * self.cl.gd.integrate(self.classical_material.charge_density)
        
        # Write 
        if self.rank == 0:
            line = '%20.8lf %20.8le %22.12le %22.12le %22.12le\n' \
                 % (self.time, norm, dm[0], dm[1], dm[2])
            self.dm_file.write(line)
            self.dm_file.flush()

    def read(self, paw, filename='poisson'):
        
        world = paw.wfs.world
        master = (world.rank == 0)
        parallel = (world.size > 1)
        self.rank = paw.wfs.world.rank

        r = gpaw_io_open(filename, 'r', world)

        version = r['version']
        parprint('reading poisson gpw-file... version: %f', version)
        self.classical_material.initialize(self.cl.gd)
        
        # Read self.classical_material.charge_density
        if self.cl.gd.comm.rank == 0:
            big_charge_density = np.array(r.get('classical_material_rho'), dtype=float)
        else:
            big_charge_density = None
        self.cl.gd.distribute(big_charge_density, self.classical_material.charge_density)
                
        # Read self.classical_material.polarization_Total
        if self.cl.gd.comm.rank == 0:
            big_polarization_total = np.array(r.get('polarization_total'), dtype=float)
        else:
            big_polarization_total = None
        self.cl.gd.distribute(big_polarization_total, self.classical_material.polarization_total)
        
        # Read self.classical_material.polarizations
        if self.cl.gd.comm.rank == 0:
            big_polarizations = np.array(r.get('polarizations'),
                                         dtype=float)
        else:
            big_polarizations = None
        self.cl.gd.distribute(big_polarizations, self.classical_material.polarizations)
        
        r.close()
        world.barrier()
         
    def write(self, paw,
                     filename='poisson'):
        # parprint('Writing FDTDPoissonSolver data to %s' % (filename))
        rho = self.classical_material.charge_density
        world = paw.wfs.world
        domain_comm = self.cl.gd.comm
        kpt_comm = paw.wfs.kpt_comm
        band_comm = paw.wfs.band_comm
    
        master = (world.rank == 0)
        parallel = (world.size > 1)
        
        w = gpaw_io_open(filename, 'w', world)
        w['history'] = 'FDTDPoissonSolver restart file'
        w['version'] = 1
        w['lengthunit'] = 'Bohr'
        w['energyunit'] = 'Hartree'
        w['DataType'] = 'Float'
        
        # Create dimensions for various netCDF variables:
        ng = self.cl.gd.get_size_of_global_array()
        
        # Write the classical charge density
        w.dimension('ngptsx', ng[0])
        w.dimension('ngptsy', ng[1])
        w.dimension('ngptsz', ng[2])
        w.add('classical_material_rho',
              ('ngptsx', 'ngptsy', 'ngptsz'),
              dtype=float,
              write=master)
        if kpt_comm.rank == 0:
            charge_density = self.cl.gd.collect(self.classical_material.charge_density)
            if master:
                w.fill(charge_density)

        # Write the total polarization
        w.dimension('3', 3)
        w.dimension('ngptsx', ng[0])
        w.dimension('ngptsy', ng[1])
        w.dimension('ngptsz', ng[2])
        w.add('polarization_total',
              ('3', 'ngptsx', 'ngptsy', 'ngptsz'),
              dtype=float,
              write=master)
        if kpt_comm.rank == 0:
            polarization_total = self.cl.gd.collect(self.classical_material.polarization_total)
            if master:
                w.fill(polarization_total)

        # Write the partial polarizations
        w.dimension('3', 3)
        w.dimension('Nj', self.classical_material.Nj)
        w.dimension('ngptsx', ng[0])
        w.dimension('ngptsy', ng[1])
        w.dimension('ngptsz', ng[2])
        w.add('polarizations',
              ('3', 'Nj', 'ngptsx', 'ngptsy', 'ngptsz'),
              dtype=float,
              write=master)
        if kpt_comm.rank == 0:
            polarizations = self.cl.gd.collect(self.classical_material.polarizations)
            if master:
                w.fill(polarizations)

        w.close()
        world.barrier()
        

