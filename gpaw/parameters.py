import numpy as np
from ase.units import Hartree, Bohr

import gpaw.mpi as mpi
from gpaw.poisson import PoissonSolver, FFTPoissonSolver
from gpaw.occupations import FermiDirac
from gpaw import parsize, parsize_bands, sl_default, sl_diagonalize, \
                 sl_inverse_cholesky, sl_lcao, buffer_size

class InputParameters(dict):
    def __init__(self, **kwargs):
        dict.__init__(self, [
            ('h',               None),  # Angstrom
            ('xc',              'LDA'),
            ('gpts',            None),
            ('kpts',            [(0, 0, 0)]),
            ('lmax',            2),
            ('charge',          0),
            ('fixmom',          False),      # don't use this
            ('nbands',          None),
            ('setups',          'paw'),
            ('basis',           {}),
            ('width',           None),  # eV, don't use this
            ('occupations',     None),
            ('spinpol',         None),
            ('usesymm',         True),
            ('stencils',        (3, 3)),
            ('fixdensity',      False),
            ('mixer',           None),
            ('txt',             '-'),
            ('hund',            False),
            ('random',          False),
            ('dtype',           float),
            ('maxiter',         120),
            ('parallel',        {'domain':              parsize,
                                 'band':                parsize_bands,
                                 'stridebands':         False,
                                 'sl_default':          sl_default,
                                 'sl_diagonalize':      sl_diagonalize,
                                 'sl_inverse_cholesky': sl_inverse_cholesky,
                                 'sl_lcao':             sl_lcao,
                                 'buffer_size':         buffer_size}),
            ('parsize',         None),
            ('parsize_bands',   None),
            ('parstride_bands', False),
            ('external',        None),  # eV
            ('verbose',         0),
            ('eigensolver',     None),
            ('poissonsolver',   None),
            ('communicator' ,   mpi.world),
            ('idiotproof'   ,   True),
            ('mode',            'fd'),
            ('convergence',     {'energy':      0.0005,  # eV / electron
                                 'density':     1.0e-4,
                                 'eigenstates': 1.0e-9,  # XXX ???
                                 'bands':       'occupied'}),
            ])
        dict.update(self, kwargs)

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        repr = 'InputParameters(**%s)' % dictrepr
        return repr
    
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        assert key in self
        self[key] = value

    def update(self, parameters):
        assert isinstance(parameters, dict)

        for key, value in parameters.items():
            assert key in self

            haschanged = (self[key] != value)

            if isinstance(haschanged, np.ndarray):
                haschanged = haschanged.any()

            #if haschanged:
            #    self.notify(key)

        dict.update(self, parameters)

    def read(self, reader):
        """Read state from file."""

        if isinstance(reader, str):
            reader = gpaw.io.open(reader, 'r')

        r = reader

        master = (reader.comm.rank == 0) # read only on root of reader.comm

        if hasattr(reader, 'hdf5'): 
            hdf5 = True
        else:
            hdf5 = False

        par_kwargs = {}
        if hdf5:
            par_kwargs.update({'parallel': False, 'read': master})

        version = r['version']
        
        assert version >= 0.3
    
        self.xc = r['XCFunctional']
        self.nbands = r.dimension('nbands')
        self.spinpol = (r.dimension('nspins') == 2)

        dim = 3 # k-point grid dimensions
        nbzkpts = r.dimension('nbzkpts')
 
        if r.has_array('NBZKPoints'):
            self.kpts = np.empty(dim, int)
            # Read on master, then broadcast
            if master:
                self.kpts = r.get('NBZKPoints', **par_kwargs)
        else:
            self.kpts = np.empty((nbzkpts, dim), float)
            # Read on master, then broadcast
            if master:
                self.kpts = r.get('BZKPoints', **par_kwargs)
        r.comm.broadcast(self.kpts, 0)

        self.usesymm = r['UseSymmetry']
        try:
            self.basis = r['BasisSet']
        except KeyError:
            pass

        if version >= '0.9':
            h = r['GridSpacing']
        else:
            h = None

        if h is None:
        ## gpts modified to account for boundary condition in non-PBC
            self.gpts = ((r.dimension('ngptsx') + 1) // 2 * 2,
                         (r.dimension('ngptsy') + 1) // 2 * 2,
                         (r.dimension('ngptsz') + 1) // 2 * 2)
        else:
            self.h = Bohr * h

        self.lmax = r['MaximumAngularMomentum']
        self.setups = r['SetupTypes']
        self.fixdensity = r['FixDensity']
        if version <= 0.4:
            # Old version: XXX
            print('# Warning: Reading old version 0.3/0.4 restart files ' +
                  'will be disabled some day in the future!')
            self.convergence['eigenstates'] = r['Tolerance']
        else:
            nbtc = r['NumberOfBandsToConverge']
            if not isinstance(nbtc, (int, str)):
                # The string 'all' was eval'ed to the all() function!
                nbtc = 'all'
            self.convergence = {'density': r['DensityConvergenceCriterion'],
                                'energy':
                                r['EnergyConvergenceCriterion'] * Hartree,
                                'eigenstates':
                                r['EigenstatesConvergenceCriterion'],
                                'bands': nbtc}
            if version <= 0.6:
                mixer = 'Mixer'
                weight = r['MixMetric']
            elif version <= 0.7:
                mixer = r['MixClass']
                weight = r['MixWeight']
                metric = r['MixMetric']
                if metric is None:
                    weight = 1.0
            else:
                mixer = r['MixClass']
                weight = r['MixWeight']

            if mixer == 'Mixer':
                from gpaw.mixer import Mixer
            elif mixer == 'MixerSum':
                from gpaw.mixer import MixerSum as Mixer
            elif mixer == 'MixerSum2':
                from gpaw.mixer import MixerSum2 as Mixer
            elif mixer == 'MixerDif':
                from gpaw.mixer import MixerDif as Mixer
            elif mixer == 'DummyMixer':
                from gpaw.mixer import DummyMixer as Mixer
            else:
                Mixer = None

            if Mixer is None:
                self.mixer = None
            else:
                self.mixer = Mixer(r['MixBeta'], r['MixOld'], weight)
            
        if version == 0.3:
            # Old version: XXX
            print('# Warning: Reading old version 0.3 restart files is ' +
                  'dangerous and will be disabled some day in the future!')
            self.stencils = (2, 3)
            self.charge = 0.0
            fixmom = False
        else:
            self.stencils = (r['KohnShamStencil'],
                             r['InterpolationStencil'])
            if r['PoissonStencil'] == 999:
                self.poissonsolver = FFTPoissonSolver()
            else:
                self.poissonsolver = PoissonSolver(nn=r['PoissonStencil'])
            self.charge = r['Charge']
            fixmom = r['FixMagneticMoment']

        self.occupations = FermiDirac(r['FermiWidth'] * Hartree,
                                      fixmagmom=fixmom)

        try:
            self.mode = r['Mode']
        except KeyError:
            self.mode = 'fd'

        try:
            dtype = r['DataType']
            if dtype=='Float':
                self.dtype = float
            else:
                self.dtype = complex
        except KeyError:
            self.dtype = float
