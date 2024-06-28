from ase.constraints import ExpCellFilter
from ase.optimize import BFGS
from gpaw import GPAW


def relax(atoms, calculator):
    atoms.calc = GPAW(**calculator)
    opt = BFGS(ExpCellFilter(atoms), trajectory='opt.traj',
               logfile='opt.log')
    opt.run(fmax=0.01)
    # Remove the calculator before returning the atoms,
    # because the calculator object as such cannot be saved:
    atoms.calc = None
    return atoms

# --- end-snippet-1 ---

def groundstate(atoms, calculator):
    from pathlib import Path
    atoms.calc = GPAW(**calculator)
    atoms.get_potential_energy()
    path = Path('groundstate.gpw')
    atoms.calc.write(path)
    return path

# --- literalinclude-divider-2 ---

def bandstructure(gpw):
    gscalc = GPAW(gpw)
    atoms = gscalc.get_atoms()
    bandpath = atoms.cell.bandpath(npoints=100)
    bandpath.write('bandpath.json')
    calc = gscalc.fixed_density(
        kpts=bandpath.kpts, symmetry='off', txt='bs.txt')
    bs = calc.band_structure()
    bs.write('bs.json')
    return bs
