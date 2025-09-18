from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[48].position + (1.0, 0.0, 1.0)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.9, 0.0, -0.4)])
write('restart.json', atoms)