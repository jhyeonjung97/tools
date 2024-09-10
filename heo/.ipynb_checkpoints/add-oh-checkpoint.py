from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[17].position + (0.0, 0.0, 2.0)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
write('restart-oh.json', atoms)