from ase.io import read, write

atoms = read('restart.json')
mag_moments = [0] * len(atoms)
mag_moments[0] = mag_moments[1] = mag_moments[2] = mag_moments[3] = 2
mag_moments[4] = mag_moments[5] = mag_moments[6] = mag_moments[7] = -2
atoms.set_initial_magnetic_moments(mag_moments)

write('restart-mag.json', atoms)