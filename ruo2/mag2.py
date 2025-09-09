from ase.io import read, write

atoms = read('restart.json')

mag_moments = [0] * len(atoms)

for i in [0, 1, 6, 7, 8, 9, 14, 15]:
    mag_moments[i] = -2
for i in [2, 3, 4, 5, 10, 11, 12, 13]:
    mag_moments[i] = 2
atoms.set_initial_magnetic_moments(mag_moments)

write('restart.json', atoms)