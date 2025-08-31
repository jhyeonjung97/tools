import numpy as np
from ase.io import read, write
from ase.build import surface
from ase.constraints import FixAtoms

atoms = read('restart.json')

atoms = surface(atoms, indices=(1,0,0), layers=3, vacuum=1.0)


# Delete lowest 10 and highest 14 atoms by z-coordinate
z_positions = atoms.get_positions()[:, 2]
order = np.argsort(z_positions)
lowest_indices = order[:8]
highest_indices = order[-16:]
indices_to_delete = np.concatenate((lowest_indices, highest_indices))
for idx in sorted(set(indices_to_delete.tolist()), reverse=True):
    del atoms[idx]

x_coords = atoms.positions[:, 0]
y_coords = atoms.positions[:, 1]
z_coords = atoms.positions[:, 2]
# First sort by z, then y, then x
sorted_indices = np.lexsort((x_coords, y_coords, z_coords))
atoms = atoms[sorted_indices]

# Fix lowest 24 atoms by z-coordinate and sort by element order Re, Ru, O
z_positions = atoms.get_positions()[:, 2]
num_to_fix = min(24, len(atoms))
fix_mask = np.zeros(len(atoms), dtype=bool)
fix_mask[np.argsort(z_positions)[:num_to_fix]] = True

# Sort atoms by desired element order
order_map = {'Re': 0, 'Ru': 1, 'O': 2}
symbols = atoms.get_chemical_symbols()
sorted_indices = sorted(range(len(atoms)), key=lambda i: order_map.get(symbols[i], 999))
atoms = atoms[sorted_indices]

# Recompute fixed indices after reordering
new_fix_mask = fix_mask[sorted_indices]
fix_indices = np.nonzero(new_fix_mask)[0].tolist()
atoms.set_constraint(FixAtoms(indices=fix_indices))

atoms.center(vacuum=9.0, axis=2)

write('slab.json', atoms)