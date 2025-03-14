from ase.io import read, write

# Read atomic structure from CONTCAR
atoms = read('CONTCAR')

# Modify z-coordinate for N and Co atoms
for atom in atoms:
    if atom.symbol == 'N':
        atom.position[2] = 10.0  # Modify z-coordinate
    elif atom.symbol == 'Co':
        dz = 10.0 - atom.position[2]
        atom.position[2] = 10.0  # Modify z-coordinate
    elif atom.symbol in ['O', 'H']:
        atom.position[2] += dz
        
# Write the modified structure to POSCAR
write('POSCAR', atoms)