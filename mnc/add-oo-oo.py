from ase.io import read, write
from ase import Atoms

atoms = read('CONTCAR')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('O', positions=[atoms[-1].position + (0.0, 0.0, -2.0)])
    atoms += Atoms('O', positions=[atoms[-1].position + (0.0, 1.0, -1.0)])
    atoms += Atoms('O', positions=[atoms[-3].position + (0.0, 0.0, 2.0)])
    atoms += Atoms('O', positions=[atoms[-1].position + (-0.1, -1.3, 0.7)])
    write('start-oo-oo.traj', atoms)