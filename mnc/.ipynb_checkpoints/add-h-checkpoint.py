from ase.io import read, write
from ase import Atoms

atoms = read('POSCAR')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 0.0, 2.0)])
    write('start.traj', atoms)