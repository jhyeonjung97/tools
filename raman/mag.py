from ase.io import read, write

atoms = read('restart.json')

initial_magnetic_moments = atoms.get_initial_magnetic_moments()
for atom in atoms:
    if atom.symbol == 'Co':
        initial_magnetic_moments[atom.index] = 2
    elif atom.symbol == 'Fe':
        initial_magnetic_moments[atom.index] = 3
    else:
        initial_magnetic_moments[atom.index] = 0
atoms.set_initial_magnetic_moments(initial_magnetic_moments)

write('restart.json', atoms)