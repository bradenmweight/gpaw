class Contribution:
    def __init__(self, weight):
        self.weight = weight

    def get_name(self):
        raise NotImplementedError

    def get_desc(self):
        raise NotImplementedError

    def initialize(self, density, hamiltonian, wfs):
        self.density = density
        self.gd = density.gd
        self.finegd = density.finegd
        self.nspins = wfs.nspins
        self.wfs = wfs

    def initialize_1d(self, ae):
        self.ae = ae

    def initialize_from_atomic_orbitals(self, basis_functions):
        # Pass if contribution needs only density which is already initialized
        pass

    def calculate_spinpaired(self, e_g, n_g, v_g):
        raise NotImplementedError

    def calculate_spinpolarized(self, e_g, na_g, va_g, nb_g, vb_g):
        raise NotImplementedError

    def calculate_energy_and_derivatives(self, setup, D_sp, H_sp):
        raise NotImplementedError

    def add_smooth_xc_potential_and_energy_1d(self, vt_g):
        raise NotImplementedError

    def add_extra_setup_data(self, dict):
        pass

    def write(self, writer):
        pass

    def read(self, reader):
        pass
