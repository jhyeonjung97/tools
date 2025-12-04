from ase.io import read, write
from ase import Atoms

atoms = read('restart.json')
atoms += Atoms('O', positions=[atoms[4].position - (0.0, 0.0, 2.0)])
atoms += Atoms('O', positions=[atoms[5].position - (0.0, 0.0, 2.0)])
atoms += Atoms('O', positions=[atoms[16].position + (0.0, 0.0, 2.0)])
atoms += Atoms('O', positions=[atoms[17].position + (0.0, 0.0, 2.0)])

min_z = atoms.positions[:,2].min()
max_z = atoms.positions[:,2].max()
height = max_z - min_z + 18.0

l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.cell = (l1, l2, height, a1, a2, a3)

atoms.center(axis=2)

atoms.write('o.json')