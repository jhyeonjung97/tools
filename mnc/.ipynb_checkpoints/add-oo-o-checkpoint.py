from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('O', positions=[atoms[-1].position + (0.1, 0.2, 2.2)])
    atoms += Atoms('O', positions=[atoms[-1].position + (-0.1, -1.3, 0.7)])
    atoms += Atoms('O', positions=[atoms[-3].position + (0.0, 0.0, -2.0)])
    write('restart-oo-o.json', atoms)