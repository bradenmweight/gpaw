import numpy as np
from ase.units import Hartree
from Scientific.IO.NetCDF import NetCDFFile
from ase.lattice.spacegroup import Spacegroup

from gpaw.wavefunctions.pw import PWDescriptor


class ETSFWriter:
    def __init__(self, filename='gpaw', title='gpaw'):
        if not filename.endswith('-etsf.nc'):
            if filename.endswith('.nc'):
                filename = filename[:-3] + '-etsf.nc'
            else:
                filename = filename + '-etsf.nc'
            
        self.nc = NetCDFFile(filename, 'w')

        self.nc.file_format = 'ETSF Nanoquanta'
        self.nc.file_format_version = np.array([3.3], dtype=np.float32)
        self.nc.Conventions = 'http://www.etsf.eu/fileformats/'
        self.nc.history = 'File generated by GPAW'
        self.nc.title = title

    def write(self, calc, spacegroup=1):

        #sg = Spacegroup(spacegroup)
        #print sg
        
        wfs = calc.wfs
        setups = wfs.setups
        bd = wfs.bd
        kd = wfs.kd
        
        atoms = calc.atoms
        natoms = len(atoms)
        
        if wfs.symmetry is None:
            op_scc = np.eye(3, dtype=int).reshape((1, 3, 3))
        else:
            op_scc = wfs.symmetry.op_scc

        specie_a = np.empty(natoms, np.int32)
        nspecies = 0
        species = {}
        names = []
        symbols = []
        numbers = []
        charges = []
        for a, id in enumerate(setups.id_a):
            if id not in species:
                species[id] = nspecies
                nspecies += 1
                names.append(setups[a].symbol)
                symbols.append(setups[a].symbol)
                numbers.append(setups[a].Z)
                charges.append(setups[a].Nv)
            specie_a[a] = species[id]
            
        dimensions = [
            ('character_string_length', 80),
            ('max_number_of_states', bd.nbands),
            ('number_of_atoms', len(atoms)),
            ('number_of_atom_species', nspecies),
            ('number_of_cartesian_directions', 3),
            ('number_of_components', 1),
            ('number_of_kpoints', kd.nibzkpts),
            ('number_of_reduced_dimensions', 3),
            ('number_of_spinor_components', 1),
            ('number_of_spins', wfs.nspins),
            ('number_of_symmetry_operations', len(op_scc)),
            ('number_of_vectors', 3),
            ('real_or_complex_coefficients', 2),
            ('symbol_length', 2)]

        for name, size in dimensions:
            print('%-34s %d' % (name, size))
            self.nc.createDimension(name, size)

        var = self.add_variable
        
        var('space_group', (), np.array(spacegroup, dtype=int))
        var('primitive_vectors',
            ('number_of_vectors', 'number_of_cartesian_directions'),
            wfs.gd.cell_cv, units='atomic units')
        var('reduced_symmetry_matrices',
            ('number_of_symmetry_operations',
             'number_of_reduced_dimensions', 'number_of_reduced_dimensions'),
            op_scc.astype(np.int32), symmorphic='yes')
        var('reduced_symmetry_translations',
            ('number_of_symmetry_operations', 'number_of_reduced_dimensions'),
            np.zeros((len(op_scc), 3), dtype=np.int32))
        var('atom_species', ('number_of_atoms',), specie_a + 1)
        var('reduced_atom_positions',
            ('number_of_atoms', 'number_of_reduced_dimensions'),
            atoms.get_scaled_positions())
        var('atomic_numbers', ('number_of_atom_species',),
            np.array(numbers, dtype=float))
        var('valence_charges', ('number_of_atom_species',),
            np.array(charges, dtype=float))
        var('atom_species_names',
            ('number_of_atom_species', 'character_string_length'), names)
        var('chemical_symbols', ('number_of_atom_species', 'symbol_length'),
            symbols)
        var('pseudopotential_types',
            ('number_of_atom_species', 'character_string_length'),
            ['HGH'] * nspecies)
        var('fermi_energy', (), calc.occupations.fermilevel,
            units='atomic units')
        var('smearing_scheme', ('character_string_length',), 'fermi-dirac')
        var('smearing_width', (), calc.occupations.width, units='atomic units')
        var('number_of_states', ('number_of_spins', 'number_of_kpoints'),
            np.zeros((wfs.nspins, kd.nibzkpts), np.int32) + bd.nbands,
            k_dependent='no')
        var('eigenvalues',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'),
            np.array([[calc.get_eigenvalues(k, s) / Hartree
                       for k in range(kd.nibzkpts)]
                      for s in range(wfs.nspins)]), units='atomic units')
        var('occupations',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'),
            np.array([[calc.get_occupation_numbers(k, s) / kd.weight_k[k]
                       for k in range(kd.nibzkpts)]
                      for s in range(wfs.nspins)]))
        var('reduced_coordinates_of_kpoints',
            ('number_of_kpoints', 'number_of_reduced_dimensions'), kd.ibzk_kc)
        var('kpoint_weights', ('number_of_kpoints',), kd.weight_k)
        var('basis_set', ('character_string_length',), 'plane_waves')
        var('kinetic_energy_cutoff', (), 1.0 * ecut, units='atomic units')
        var('number_of_coefficients', ('number_of_kpoints',),
            np.zeros(kd.nibzkpts, np.int32) + len(i_Gc),
            k_dependent='no')
        var('number_of_electrons', (), np.array(wfs.nvalence, dtype=np.int32))
        self.nc.close()

    def add_variable(self, name, dims, data=None, **kwargs):
        if data is None:
            char = 'd'
        else:
            if isinstance(data, np.ndarray):
                char = data.dtype.char
            elif isinstance(data, float):
                char = 'd'
            elif isinstance(data, int):
                char = 'i'
            else:
                char = 'c'
        print('%-34s %s%s' % (
            name, char,
            tuple([self.nc.dimensions[dim] for dim in dims])))
        var = self.nc.createVariable(name, char, dims)
        for attr, value in kwargs.items():
            setattr(var, attr, value)
        if data is not None:
            if len(dims) == 0:
                var.assignValue(data)
            else:
                if char == 'c':
                    if len(dims) == 1:
                        var[:len(data)] = data
                    else:
                        for i, x in enumerate(data):
                            var[i, :len(x)] = x
                else:
                    var[:] = data
        return var
