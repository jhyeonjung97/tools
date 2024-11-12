from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('H', positions=[atoms[-2].position + (0.0, 0.0, 1.2)])
    write('restart-nh.json', atoms)