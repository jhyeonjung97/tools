from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[48].position + (1.0, 0.0, 1.0)])
atoms += Atoms('O', positions=[atoms[49].position + (-0.6, 0.0, 0.8)])
write('restart.json', atoms)