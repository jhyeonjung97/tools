from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('O', positions=[atoms[-1].position + (0.2, -0.1, 2.2)])
    atoms += Atoms('O', positions=[atoms[-1].position + (-1.3, 0.1, 0.7)])
    atoms += Atoms('H', positions=[atoms[-1].position + (0.1, -0.8, 0.5)])
    write('restart-ooh.json', atoms)