from ase.io import read, write

atoms = read('restart.json')
z_coords = atoms.get_positions()[:, 2]
indices_to_delete = z_coords.argsort()[9]

for index in indices_to_delete:
    del atoms[index]
write('restart_del.json', atoms)