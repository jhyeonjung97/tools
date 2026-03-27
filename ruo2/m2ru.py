from ase.io import read, write
atoms = read('restart.json')

atoms[19].symbol = 'Ru'

magmoms = atoms.get_initial_magnetic_moments()
magmoms[19] = magmoms[16]

atoms.set_initial_magnetic_moments(magmoms)

write('restart.json', atoms)