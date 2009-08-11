from ase import Hartree
import numpy as np

def collect_D_asp(density):
    all_D_asp = []
    for a, setup in enumerate(density.setups):
        D_sp = density.D_asp.get(a)
        if D_sp is None:
            ni = setup.ni
            D_sp = np.empty((density.nspins, ni * (ni + 1) // 2))
        if density.gd.comm.size > 1:
            density.gd.comm.broadcast(D_sp, density.rank_a[a])
        all_D_asp.append(D_sp)      
    return all_D_asp

def distribute_D_asp(D_asp, density):
    for a in range(len(density.setups)):
        if density.D_asp.get(a) is not None:
            density.D_asp[a] = D_asp[a]
   
class Side:
    def __init__(self, type, atoms, direction):
        self.type = type
        self.atoms = atoms
        self.direction = direction
        self.n_atoms = len(atoms)
        calc = atoms.calc
        self.N_c = calc.gd.N_c

    def abstract_boundary(self):
        calc = self.atoms.calc
        gd = calc.gd
        finegd = calc.finegd
        nn = finegd.N_c[2]
        ns = calc.wfs.nspins

        dim = gd.N_c
        d1 = dim[0] // 2
        d2 = dim[1] // 2
        
        vHt_g = finegd.collect(calc.hamiltonian.vHt_g, True)
        self.boundary_vHt_g = self.slice(nn, vHt_g)
        
        vt_sg = finegd.collect(calc.hamiltonian.vt_sg, True)
        self.boundary_vt_sg_line = self.slice(nn, vt_sg[:, d1 * 2, d2 * 2])
        
        nt_sg = finegd.collect(calc.density.nt_sg, True)
        self.boundary_nt_sg = self.slice(nn, nt_sg)        
        
        rhot_g = finegd.collect(calc.density.rhot_g, True)
        self.boundary_rhot_g_line = self.slice(nn, rhot_g[d1 * 2, d2 * 2])
  
        nn /= 2
        vt_sG = gd.collect(calc.hamiltonian.vt_sG, True)
        self.boundary_vt_sG_line = self.slice(nn, vt_sG[:, d1, d2])
        
        nt_sG = calc.gd.collect(calc.density.nt_sG, True)
        self.boundary_nt_sG = self.slice(nn, nt_sG)
        
        self.D_asp = collect_D_asp(calc.density)
        
    def slice(self, nn, in_array):
        if self.type == 'LR':
            seq1 = np.arange(-nn + 1, 1)
            seq2 = np.arange(nn)
            di = len(in_array.shape) - 1
            if self.direction == '-':
                slice_array = np.take(in_array, seq1, axis=di)
            else:
                slice_array = np.take(in_array, seq2, axis=di)
        return slice_array

class Surrounding:
    def __init__(self, tp, type='LR'):
        self.tp = tp
        self.type = type
        self.lead_num = tp.lead_num
        self.initialize()
        
    def initialize(self):
        if self.type == 'LR':
            self.sides = {}
            self.bias_index = {}
            self.side_basis_index = {}
            self.nn = []
            self.directions = ['-', '+']
            for i in range(self.lead_num):
                direction = self.directions[i]
                side = Side('LR', self.tp.atoms_l[i], direction)
                self.sides[direction] = side
                self.bias_index[direction] = self.tp.bias[i]
                self.side_basis_index[direction] = self.tp.lead_index[i]                
                self.nn.append(side.N_c[2])
            self.nn = np.array(self.nn)
            self.operator = self.tp.hamiltonian.poisson.operators[0]
            
        elif self.type == 'all':
            raise NotImplementError()
        self.calculate_sides()
        self.initialized = True

    def reset_bias(self, bias):
        self.bias = bias
        for i in range(self.lead_num):
            direction = self.directions[i]
            self.bias_index[direction] = bias[i]
        self.combine()

    def calculate_sides(self):
        if self.type == 'LR':
            for name, in self.sides:
                self.sides[name].abstract_boundary()
        if self.type == 'all':
            raise NotImplementError('type all not yet')
            
    def get_extra_density(self):
        if self.type == 'LR':
            rhot_g = self.tp.finegd.zeros()
            self.operator.apply(self.tp.hamiltonian.vHt_g, rhot_g)
            nn = self.nn[0] * 2
            self.extra_rhot_g = self.uncapsule(nn, rhot_g)

    def capsule(self, nn, loc_in_array):
        ns = self.tp.nspins
        gd, gd0 = self.tp.finegd, self.tp.finegd0
        cap_array = gd.collect(self.tp.hamiltonian.vHt_g, True)
        in_array = gd0.collect(loc_in_array, True)
        if len(in_array.shape) == 4:
            local_cap_array = gd.zeros(ns)
            cap_array[:, :, :, nn:-nn] = in_array
        else:
            local_cap_array = gd.zeros()
            cap_array[:, :, nn:-nn] = in_array
        gd.distribute(cap_array, local_cap_array)
        return local_cap_array
    
    def uncapsule(self, nn, in_array, nn2=None):
        gd, gd0 = self.tp.finegd, self.tp.finegd0
        nn1 = nn
        if nn2 == None:
            nn2 = nn1
        di = 2
        local_uncap_array = gd0.zeros()
        global_in_array = gd.collect(in_array, True)
        seq = np.arange(nn1, global_in_array.shape[di] - nn2)    
        uncap_array = np.take(global_in_array, seq, axis=di)
        gd0.distribute(uncap_array, local_uncap_array)
        return local_uncap_array
      
    def combine(self):
        if self.type == 'LR':
            nn = self.nn[0] * 2
            ham = self.tp.hamiltonian
            if ham.vt_sg is None:
                ham.vt_sg = ham.finegd.empty(ham.nspins)
                ham.vHt_g = ham.finegd.zeros()
                ham.vt_sG = ham.gd.empty(ham.nspins)
                ham.poisson.initialize()
            vHt_g = ham.finegd.zeros(global_array=True)
            bias_shift0 = self.bias_index['-'] / Hartree
            bias_shift1 = self.bias_index['+'] / Hartree
            vHt_g[:, :, :nn] = self.sides['-'].boundary_vHt_g + bias_shift0
            vHt_g[:, :, -nn:] = self.sides['+'].boundary_vHt_g + bias_shift1
            ham.finegd.distribute(vHt_g, ham.vHt_g)
            self.get_extra_density()

    def combine_vHt_g(self, vHt_g):
        nn = self.nn[0] * 2
        self.tp.hamiltonian.vHt_g = self.capsule(nn, vHt_g)

    def combine_nt_sG(self):
        nn = self.nn[0]
        gd = self.tp.gd
        nt_sG = gd.collect(self.tp.density.nt_sG, True)
        nt_sG[:, :, :, :nn] = self.sides['-'].boundary_nt_sG
        nt_sG[:, :, :, -nn:] = self.sides['+'].boundary_nt_sG
        gd.distribute(nt_sG, self.tp.density.nt_sG)
        
    def combine_nt_sg(self):
        nn = self.nn[0] * 2
        gd = self.tp.finegd
        nt_sg = gd.collect(self.tp.density.nt_sg, True)
        nt_sg[:, :, :, :nn] = self.sides['-'].boundary_nt_sg
        nt_sg[:, :, :, -nn:] = self.sides['+'].boundary_nt_sg
        gd.distribute(nt_sg, self.tp.density.nt_sg)        

    def combine_D_asp(self):
        density = self.tp.density
        all_D_asp = collect_D_asp(density)
        nao = len(self.tp.original_atoms)
        for i in range(self.lead_num):
            direction = self.directions[i]
            side = self.sides[direction]
            for n in range(side.n_atoms):
                all_D_asp[nao + n] = side.D_asp[n]
            nao += side.n_atoms
        distribute_D_asp(all_D_asp, density)     
       
    def abstract_inner_rhot(self):
        nn = self.nn[0] * 2
        rhot_g = self.uncapsule(nn, self.tp.density.rhot_g)
        rhot_g -= self.extra_rhot_g
        return rhot_g
        