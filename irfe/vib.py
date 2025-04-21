from ase.io import read, write
from ase.constraints import FixAtoms

# Read the structure
atoms = read('CONTCAR')

# Create a mask for metal atoms (Ir, Mn, Fe, Co, Ni)
metal_symbols = ['Ir', 'Mn', 'Fe', 'Co', 'Ni']
metal_indices = [i for i, atom in enumerate(atoms) if atom.symbol in metal_symbols]

# Apply constraints to fix metal atoms
constraint = FixAtoms(indices=metal_indices)
atoms.set_constraint(constraint)

# Write the structure with constraints
write('vib.vasp', atoms)