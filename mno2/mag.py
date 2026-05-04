import os
from ase.io import read, write

for file in os.listdir('.'):
    if file.endswith('.vasp'):
        atoms = read(file)
        mag_moments = [0] * len(atoms)
        for atom in atoms:
            if atom.symbol == 'Mn':
                mag_moments[atom.index] = 3
        atoms.set_initial_magnetic_moments(mag_moments)
        write(file.replace('.vasp', '.json'), atoms)