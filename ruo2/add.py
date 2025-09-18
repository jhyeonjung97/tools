from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[48].position + (0.0, -0.8, 1.2)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.0, -0.9, -0.4)])
write('restart.json', atoms)