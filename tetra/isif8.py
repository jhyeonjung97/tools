import os
import glob
from ase.io import read, write

# Define the base directory pattern
base_path = "/pscratch/sd/j/jiuy97/3_V_bulk/isif8/*_*_*/*d/*_*/"

def cell_shape(atoms, coord, output_filename):
    a = atoms.cell.lengths()[0]
    if coord == 'WZ':
        new_lengths = [a, a, a * (2 * (6**0.5) / 3)]
        new_angles = [90, 90, 120]
    elif coord == 'ZB':
        new_lengths = [a, a, a]
        new_angles = [90, 90, 90]
    elif coord == 'TN':
        new_lengths = [a, a * (3.42 / 4.68), a * (5.13 / 4.68)]
        new_angles = [90, 99.54, 90]
    elif coord == 'PD':
        new_lengths = [a, a, a * (2 * (6**0.5) / 3)]
        new_angles = [90, 90, 90]
    elif coord == 'NB':
        new_lengths = [a, a, a]
        new_angles = [60, 60, 60]
    elif coord == 'RS':
        new_lengths = [a, a, a]
        new_angles = [90, 90, 90]
    elif coord == 'LT':
        new_lengths = [a, a, a * (2**0.5)]
        new_angles = [90, 90, 90]
    elif coord == 'AQ':
        new_lengths = [a, a * (5.49 / 5.90), a * (4.75 / 5.90)]
        new_angles = [90, 90, 90]
    elif coord == 'AU':
        new_lengths = [a, a, a * (3**0.5)]
        new_angles = [90, 90, 120]
    else:
        raise ValueError(f"Unknown coordination type: {coord}")

    # Set the new cell
    atoms.set_cell([new_lengths, new_angles])
    # Write to the specified output file
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

    print(f"Processed {filename} in {dir_path} and saved to {output_filename}")
