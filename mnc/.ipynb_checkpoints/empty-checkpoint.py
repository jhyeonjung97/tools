from ase.io import read, write
from ase import Atoms

atoms = read('./1_LS/relaxed/restart.json')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('O', positions=[atoms[-1].position + (0.0, 0.0, 2.0)])
    atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
    write('restart.json', atoms)