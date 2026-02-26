from typing import Any


from ase.io import read, write
from ase.geometry.geometry import get_duplicate_atoms

atoms = read('POSCAR')
get_duplicate_atoms(atoms, cutoff=0.1, delete=True)

symbols = atoms.get_chemical_symbols()
ba_indices = [i for i, sym in enumerate(symbols) if sym == 'Ba']
co_indices = [i for i, sym in enumerate(symbols) if sym == 'Co']
fe_indices = [i for i, sym in enumerate(symbols) if sym == 'Fe']
zr_indices = [i for i, sym in enumerate(symbols) if sym == 'Zr']
y_indices = [i for i, sym in enumerate(symbols) if sym == 'Y']
o_indices = [i for i, sym in enumerate(symbols) if sym == 'O']
h_indices = [i for i, sym in enumerate(symbols) if sym == 'H']
new_indices = ba_indices + co_indices + fe_indices + zr_indices + y_indices + o_indices + h_indices
atoms = atoms[new_indices]

# l1 = atoms.cell.lengths()[0]
# l2 = atoms.cell.lengths()[1]
# # l3 = atoms.cell.lengths()[2]
# a1 = atoms.cell.angles()[0]
# a2 = atoms.cell.angles()[1]
# a3 = atoms.cell.angles()[2]
# atoms.cell = (l1, l2, 30, a1, a2, a3)

write('start.traj',atoms)