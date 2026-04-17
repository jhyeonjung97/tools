from ase.io import read, write

atoms = read('restart.json')
z_coords = atoms.get_positions()[:, 2]
indices_to_delete = z_coords.argsort()[9]

del atoms[indices_to_delete]
write('restart_del.json', atoms)