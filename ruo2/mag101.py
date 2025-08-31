from ase.io import read, write

atoms = read('restart.json')
# Set magnetic moments: 0,1,4,5 -> 2, 2,3,6,7 -> -2, rest -> 0
mag_moments = [0] * len(atoms)
mag_moments[0] = mag_moments[1] = mag_moments[4] = mag_moments[5] = 2
mag_moments[2] = mag_moments[3] = mag_moments[6] = mag_moments[7] = -2
atoms.set_initial_magnetic_moments(mag_moments)

write('restart-mag.json', atoms)