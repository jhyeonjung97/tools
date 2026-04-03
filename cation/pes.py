from ase.io import read, write

atoms = read('restart.json')

z_coords = atoms.get_positions()[:, 2]
min_z_index = z_coords.argmin()

for i in range(10):
    atoms[min_z_index].z -= 0.1
    write(f'restart_{i+1:02d}.json', atoms)
