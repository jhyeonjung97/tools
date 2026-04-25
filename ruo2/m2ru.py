from ase.io import read, write
atoms = read('restart.json')

atoms[45].symbol = 'Ru'

magmoms = atoms.get_initial_magnetic_moments()
magmoms[45] = magmoms[44]

atoms.set_initial_magnetic_moments(magmoms)

write('restart-m2ru.json', atoms)