from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np

# Read the structure
atoms = read('CONTCAR')

# Get z-coordinates of all atoms
z_coords = atoms.get_positions()[:, 2]

# Find indices of the 24 atoms with lowest z-coordinates
bottom_24_indices = np.argsort(z_coords)[:24]

# Apply constraints to fix the bottom 24 atoms
constraint = FixAtoms(indices=bottom_24_indices)
atoms.set_constraint(constraint)

# Write the structure with constraints
write('vib.vasp', atoms)