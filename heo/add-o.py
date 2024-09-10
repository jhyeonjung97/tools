from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[17].position + (0.0, 0.0, 2.0)])
write('restart-o.json', atoms)