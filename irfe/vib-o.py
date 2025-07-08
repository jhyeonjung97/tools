from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np

# Read the structure
atoms = read('CONTCAR')

# Get z-coordinates of all atoms
z_coords = atoms.get_positions()[:, 2]

# Find indices of the 48 atoms with lowest z-coordinates
bottom_48_indices = np.argsort(z_coords)[:48]

# Apply constraints to fix the bottom 48 atoms
constraint = FixAtoms(indices=bottom_48_indices)
atoms.set_constraint(constraint)

# Write the structure with constraints
write('vib-o.vasp', atoms)