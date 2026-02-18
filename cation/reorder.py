from ase.io import read, write
from ase.geometry.geometry import get_duplicate_atoms

atoms = read('restart.json')
get_duplicate_atoms(atoms, cutoff=0.1, delete=True)

symbols = atoms.get_chemical_symbols()
pt_indices = [i for i, sym in enumerate(symbols) if sym == 'Pt']
au_indices = [i for i, sym in enumerate(symbols) if sym == 'Au']
cation_indices = [i for i, sym in enumerate(symbols) if sym == 'Li' or sym == 'Na' or sym == 'K' or sym == 'Rb' or sym == 'Cs' or sym == 'Fr']
o_indices = [i for i, sym in enumerate(symbols) if sym == 'O']
h_indices = [i for i, sym in enumerate(symbols) if sym == 'H']
other_indices = [i for i, sym in enumerate(symbols) if sym != 'Pt' and sym != 'Au' and sym != 'O' and sym != 'H' and sym != 'Li' and sym != 'Na' and sym != 'K' and sym != 'Rb' and sym != 'Cs' and sym != 'Fr']
new_indices = pt_indices + au_indices + cation_indices + o_indices + h_indices + other_indices
atoms = atoms[new_indices]

# l1 = atoms.cell.lengths()[0]
# l2 = atoms.cell.lengths()[1]
# # l3 = atoms.cell.lengths()[2]
# a1 = atoms.cell.angles()[0]
# a2 = atoms.cell.angles()[1]
# a3 = atoms.cell.angles()[2]
# atoms.cell = (l1, l2, 30, a1, a2, a3)

write('restart.json',atoms)