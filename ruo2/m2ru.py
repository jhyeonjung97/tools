from ase.io import read, write
atoms = read('restart.json')

atoms[33].symbol = 'Ru'

magmoms = atoms.get_initial_magnetic_moments()
magmoms[33] = magmoms[32]

atoms.set_initial_magnetic_moments(magmoms)

write('restart-m2ru.json', atoms)