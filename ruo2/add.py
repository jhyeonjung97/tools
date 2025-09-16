from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[48].position + (0.8, 0.0, 1.2)])
atoms += Atoms('H', positions=[atoms[49].position + (-0.6, 0.0, 0.8)])
write('restart.json', atoms)