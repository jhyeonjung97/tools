import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

atoms = read('restart.json')
z_positions = atoms.get_positions()[:, 2]
num_to_fix = min(48, len(atoms))
fix_mask = np.zeros(len(atoms), dtype=bool)
fix_mask[np.argsort(z_positions)[:num_to_fix]] = True

fix_indices = np.nonzero(fix_mask)[0].tolist()
atoms.set_constraint(FixAtoms(indices=fix_indices))
write('restart.json', atoms)