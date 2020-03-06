# -*- coding: utf-8 -*-

import os.path
import numpy as np
from ase import parallel as mpi
from ase.parallel import parprint


class FiniteDifference:
    def __init__(self, atoms, propertyfunction,
                 save=False, name='fd', ending='',
                 d=0.001, parallel=1, world=None):
        """
    atoms: Atoms object
        The atoms to work on.
    propertyfunction: function that returns a single number.
        The finite difference calculation is progressed on this value.
        For proper parallel usage the function should either be
        either a property of the atom object
            fd = FiniteDifference(atoms, atoms.property_xyz)
        or an arbitrary function with the keyword "atoms"
            fd = FiniteDifference(atoms, function_xyz)
            xyz = fd.run(atoms=atoms)
    d: float
        Magnitude of displacements.
    save: If true the write statement of the calculator is called
        to save the displacementsteps.
    name: string
        Name for restart data
    ending: string
        File handel for restart data
    parallel: int
        splits the mpi.world into 'parallel' subprocs that calculate
        displacements of different atoms individually.
    """

        self.atoms = atoms
        self.indices = np.asarray(range(len(atoms)))
        self.propertyfunction = propertyfunction
        self.save = save
        self.name = name
        self.ending = ending
        self.d = d

        self.set_parallel(parallel, world)

    def calculate(self, a, i, filename='fd', **kwargs):
        """Evaluate finite difference  along i'th axis on a'th atom.
        This will trigger two calls to propertyfunction(), with atom a moved
        plus/minus d in the i'th axial direction, respectively.
        if save is True the +- states are saved after
        the calculation
        """
        if 'atoms' in kwargs:
            kwargs['atoms'] = self.atoms

        p0 = self.atoms.positions[a, i]

        self.atoms.positions[a, i] += self.d
        eplus = self.propertyfunction(**kwargs)
        if self.save is True:
            savecalc = self.atoms.get_calculator()
            savecalc.write(filename + '+' + self.ending)

        self.atoms.positions[a, i] -= 2 * self.d
        eminus = self.propertyfunction(**kwargs)
        if self.save is True:
            savecalc = self.atoms.get_calculator()
            savecalc.write(filename + '-' + self.ending)
        self.atoms.positions[a, i] = p0

        self.value[a, i] = (eminus - eplus) / (2 * self.d)
        
        if self.parallel > 1 and self.comm.rank == 0:
            print('# rank', self.world.rank, 'Atom', a,
                  'direction', i, 'FD: ', self.value[a, i])
        else:
            parprint('Atom', a, 'direction', i,
                     'FD: ', self.value[a, i])

    def run(self, **kwargs):
        """Evaluate finite differences for all atoms
        """
        self.value = np.zeros([len(self.atoms), 3])
        
        for filename, a, i in self.displacements():
            if a in self.myindices:
                self.calculate(a, i, filename=filename, **kwargs)

        self.world.barrier()
        self.value /= self.cores_per_atom
        self.world.sum(self.value)
        
        return self.value

    def displacements(self):
        for a in self.indices:
            for i in range(3):
                filename = ('{0}_{1}_{2}'.format(self.name, a, 'xyz'[i]))
                yield filename, a, i

    def restart(self, restartfunction, **kwargs):
        """Uses restartfunction to recalculate values
        from the saved files.
        If a file with the corresponding name is found the
        restartfunction is called to get the FD value
        The restartfunction should take a string as input
        parameter like the standart read() function.
        If no file is found, a calculation is initiated.
        Example:
            def re(self, name):
                calc = Calculator(restart=name)
                return calc.get_potential_energy()

            fd = FiniteDifference(atoms, atoms.get_potential_energy)
            fd.restart(re)
        """
        for filename, a, i in self.displacements():

            if (os.path.isfile(filename + '+' + self.ending) and
                    os.path.isfile(filename + '-' + self.ending)):
                eplus = restartfunction(
                    self, filename + '+' + self.ending, **kwargs)
                eminus = restartfunction(
                    self, filename + '-' + self.ending, **kwargs)
                self.value[a, i] = (eminus - eplus) / (2 * self.d)
            else:
                self.calculate(a, i, filename=filename, **kwargs)

        return self.value

    def set_parallel(self, parallel, world):
        """Copy object onto different communicators"""
        if world is None:
            world = mpi.world
        self.world = world
        
        if parallel > world.size:
            parprint('#', (self.__class__.__name__ + ':'),
                     'Serial calculation, keyword parallel ignored.')
            parallel = 1
        self.parallel = parallel

        assert world.size % parallel == 0

        # number of atoms to calculate
        natoms = len(self.atoms)
        myn = -(-natoms // parallel)  # ceil divide

        # my workers index
        self.cores_per_atom = world.size // parallel
        myi = world.rank // self.cores_per_atom
        self.myindices = list(range(min(myi * myn, natoms),
                                    min((myi + 1) * myn, natoms)))
        
        if parallel < 2:  # no redistribution needed
            return
        
        calc = self.atoms.get_calculator()
        calc.write(self.name + '_eq' + self.ending)

        ranks = np.array(range(world.size), dtype=int)
        self.ranks = ranks.reshape(
            parallel, world.size // parallel)

        for i in range(parallel):
            if world.rank in self.ranks[i]:
                comm = world.new_communicator(self.ranks[i])
                calc2 = calc.__class__.read(
                    self.name + '_eq' + self.ending,
                    communicator=comm)
                self.atoms.set_calculator(calc2)
                self.comm = comm
