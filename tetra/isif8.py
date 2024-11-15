import os
import glob
from ase.io import read, write

# Define the base directory pattern
base_path = "/pscratch/sd/j/jiuy97/3_V_bulk/isif8/*_*_*/fm/*_*/"

def cell_shape(atoms, coord, output_filename):
    a = atoms.cell.lengths()[0]  # Base lattice parameter
    if coord == 'WZ':
        new_cell = [a, a, a * (2 * (6**0.5) / 3), 90, 90, 120]
    elif coord == 'ZB':
        new_cell = [a, a, a, 33.56, 33.56, 33.56]
    elif coord == 'TN':
        new_cell = [a, a * (3.42 / 4.68), a * (5.13 / 4.68), 90, 99.54, 90]
    elif coord == 'PD':
        new_cell = [a, a, a * (3**0.5), 90, 90, 90]
    elif coord == 'NB':
        new_cell = [a, a, a, 60, 60, 60]
    elif coord == 'RS':
        new_cell = [a, a, a, 33.56, 33.56, 33.56]
    elif coord == 'LT':
        new_cell = [a, a, a * (2**0.5), 90, 90, 90]
    elif coord == 'AQ':
        new_cell = [a, a * (5.49 / 5.90), a * (4.75 / 5.90), 90, 90, 90]
    elif coord == 'AU':
        new_cell = [a, a, a * (3**0.5), 90, 90, 120]
    else:
        raise ValueError(f"Unknown coordination type: {coord}")

    fractional_positions = atoms.get_scaled_positions()
    atoms.set_cell(new_cell)
    atoms.set_scaled_positions(fractional_positions)
    write(output_filename, atoms)

# Iterate through directories matching the pattern
for dir_path in glob.glob(base_path):
    # Split the directory path by '/'
    path_parts = dir_path.strip('/').split('/')
    
    # Extract specific parts of the path
    path1 = path_parts[-1]  # Last part of the path
    path2 = path_parts[-2]  # Second last part of the path
    path3 = path_parts[-3]  # Third last part of the path

    # Extract additional information from path1 and path3
    numb = path1.split('_')[0][:2]  # First two characters from the first part of path1
    tag = path1.split('_')[0][2:3]  # Third character from the first part of path1
    metal = path1.split('_')[1]  # Second part of path1
    coord = path3.split('_')[2]  # Third part of path3

    # Change to the current directory
    os.chdir(dir_path)

    # Read the atomic structure file
    if os.path.exists('./restart.json'):
        filename = 'restart.json'
    elif os.path.exists('./start.traj'):
        filename = 'start.traj'
    else:
        print(f"Error: Neither 'restart.json' nor 'start.traj' found in {dir_path}")
        continue

    # Adjust the cell shape based on the coordination
    atoms = read(filename)
    output_filename = f"updated_{filename}"  # Save with a new name
    cell_shape(atoms, coord, output_filename)