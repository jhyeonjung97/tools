from ase.io import read, write

atoms = read('xtbopt.xyz')

atoms.cell = (30, 30, 40, 90, 90, 90)
atoms.center()
write('xtbopt.json', atoms)