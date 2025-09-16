from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[14].position + (0.0, 0.0, 2.0)])
atoms += Atoms('O', positions=[atoms[15].position + (0.0, 0.0, 2.0)])
write('restart.json', atoms)