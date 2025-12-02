from ase.io import read, write

# Read atomic structure from CONTCAR
atoms = read('POSCAR')

n = 0
# Modify z-coordinate for N and Co atoms
for atom in atoms:
    if atom.symbol == 'N':
        atom.position[2] = 10.0  # Modify z-coordinate
    elif atom.symbol == 'Co':
        dz = 10.0 + n * 0.2 - atom.position[2]
        atom.position[2] = 10.0 + n * 0.2  # Modify z-coordinate
    elif atom.symbol in ['O', 'H']:
        atom.position[2] += dz
        
# Write the modified structure to POSCAR
write('start.traj', atoms)