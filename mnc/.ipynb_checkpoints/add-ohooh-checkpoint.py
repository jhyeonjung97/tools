from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    atoms += Atoms('O', positions=[atoms[-1].position + (0.0, -1.2, 2.0)])
    atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 0.6, 0.8)])
    atoms += Atoms('O', positions=[atoms[-3].position + (0.0, 1.2, 2.0)])
    atoms += Atoms('O', positions=[atoms[-1].position + (0.0, 0.0, 1.4)])
    atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 1.0, 0.0)])
    write('restart-ohooh.json', atoms)