from ase import Atoms
from gpaw import GPAW
from gpaw.mpi import world, serial_comm

a = Atoms('H2',
          [(0, 0, 0), (0, 0, 0.74)],
          cell=(3, 3, 3),
          pbc=1)

a.calc = GPAW(mode='pw',
              nbands=8,
              parallel={'band': world.size},
              dtype=complex,
              basis='dzp', txt=None)

a.get_potential_energy()
a.calc.write('H2')

if world.size == 1:
    scalapack = None
else:
    mb = {2: 1, 4: 2, 8: 2}[world.size]
    scalapack = (mb, mb, 32)

w1 = a.calc.get_pseudo_wave_function(0)
e1 = a.calc.get_eigenvalues()

a.calc.diagonalize_full_hamiltonian(nbands=120, scalapack=scalapack)
w2 = a.calc.get_pseudo_wave_function(0)
e2 = a.calc.get_eigenvalues()

calc = GPAW('H2', txt=None, parallel={'band': world.size})
calc.diagonalize_full_hamiltonian(nbands=120, scalapack=scalapack)
w3 = calc.get_pseudo_wave_function(0)
e3 = calc.get_eigenvalues()

calc.write('H2wf', 'all')

calc = GPAW('H2wf', txt=None, communicator=serial_comm)
w4 = calc.get_pseudo_wave_function(0)
e4 = calc.get_eigenvalues()

for w in [w2, w3, w4]:
    assert abs(abs(w[1, 2, 3]) - abs(w[1, 2, 3])) < 1e-10

for e in [e2, e3, e4]:
    assert abs(e[1] - e1[1]) < 1e-10
    assert abs(e[-1] - e2[-1]) < 1e-10
