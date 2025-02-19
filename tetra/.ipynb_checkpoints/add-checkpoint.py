import os
import glob
import numpy as np
from ase import Atoms
from ase.io import read, write

def add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, z_addition, o_file, oh_file):
    l1 = atoms.cell.lengths()[0]
    l2 = atoms.cell.lengths()[1]
    l3 = atoms.cell.lengths()[2]
    a1 = atoms.cell.angles()[0]
    a2 = atoms.cell.angles()[1]
    a3 = atoms.cell.angles()[2]
    atoms.cell = (l1, l2, l3+z_addition, a1, a2, a3)
    
    # Add oxygen atoms
    for position in oxygen_positions:
        atoms += Atoms('O', positions=[position])
    write(o_file, atoms)
    
    # Add hydrogen atoms
    for position in hydrogen_positions:
        atoms += Atoms('H', positions=[position])
    write(oh_file, atoms)

# Loop through directories
for dir in glob.glob('/pscratch/sd/j/jiuy97/4_V_slab/*_*_*/*d/*_*/'):
    # Skip directories with specific patterns in their names
    if any(x in dir for x in ['x_', 's_', 'z_', '8_Tetrahedral_AQ', '9_SquarePlanar_AU']):
        continue
    
    # Change to the current directory
    os.chdir(dir)

    # Check if the 'DONE' file exists
    if not os.path.exists('DONE'):
        print(f"'DONE' file missing in directory: {dir}")
        continue

    # Read the atoms object from the file
    try:
        atoms = read('restart.json')
    except Exception as e:
        print(f"Error reading 'restart.json' in {dir}: {e}")
        continue
    
    metal_indices = [atom.index for atom in atoms if atom.symbol != 'O']
    oxygen_indices = [atom.index for atom in atoms if atom.symbol == 'O']
    
    try:
        if '1_Tetrahedral_WZ' in dir:
            bond_vector = atoms[oxygen_indices[-2]].position - atoms[metal_indices[-3]].position
            oxygen_positions = [atoms[metal_indices[-1]].position + bond_vector]
            hydrogen_positions = [oxygen_positions[0] + (0.0, 1.0, 0.0)]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o.json', 'restart-oh.json')
        
        elif '2_Tetrahedral_ZB' in dir:
            bond_vector = atoms[oxygen_indices[-2]].position - atoms[metal_indices[-4]].position
            oxygen_positions = [atoms[metal_indices[-2]].position + bond_vector,
                                atoms[metal_indices[-1]].position + bond_vector]
            hydrogen_positions = [pos + (-0.7, 0.7, 0.0) for pos in oxygen_positions]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o.json', 'restart-oh.json')

        # elif '3_SquarePlanar_TN' in dir:
        #     bond_vector = atoms[oxygen_indices[-4]].position - atoms[metal_indices[-8]].position
        #     oxygen_positions = [atoms[metal_indices[-4]].position + bond_vector,
        #                         atoms[metal_indices[-3]].position + bond_vector]
        #     hydrogen_positions = [pos + (-1.0, 0.0, 0.0) for pos in oxygen_positions]
        #     add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o1.json', 'restart-oh1.json')
            
        #     # Reset atoms and perform the second operation
        #     atoms = read('restart.json')
        #     bond_vector = (0.0, 0.5, 2.0)
        #     oxygen_positions = [atoms[metal_indices[-2]].position + bond_vector,
        #                         atoms[metal_indices[-1]].position + bond_vector]
        #     hydrogen_positions = [pos + (-1.0, 0.0, 0.0) for pos in oxygen_positions]
        #     add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o2.json', 'restart-oh2.json')

        elif '4_SquarePlanar_PD' in dir:
            bond_vector = (0.7, 0.0, 2.2)
            oxygen_positions = [atoms[metal_indices[-1]].position + bond_vector]
            hydrogen_positions = [oxygen_positions[0] + (0.0, 1.0, 0.0)]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o1.json', 'restart-oh1.json')

            atoms = read('restart.json')
            if '02_Hf' in dir or '03_Ta' in dir:
                bond_vector = atoms[oxygen_indices[-2]].position - atoms[metal_indices[-3]].position
            else:      
                bond_vector = atoms[oxygen_indices[-2]].position - atoms[metal_indices[-4]].position
            oxygen_positions = [atoms[metal_indices[-2]].position + bond_vector]
            hydrogen_positions = [oxygen_positions[0] + (0.0, 1.0, 0.0)]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o2.json', 'restart-oh2.json')

        elif '5_SquarePlanar_NB' in dir:
            bond_vector = (0.0, 0.0, atoms[oxygen_indices[-1]].position[2] - atoms[metal_indices[-3]].position[2])
            oxygen_positions = [atoms[metal_indices[-2]].position + bond_vector,
                                atoms[metal_indices[-1]].position + bond_vector]
            hydrogen_positions = [pos + (-0.7, 0.7, 0.0) for pos in oxygen_positions]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o.json', 'restart-oh.json')

        elif '6_Octahedral_RS' in dir:
            bond_vector = (0.0, 0.0, atoms[oxygen_indices[-1]].position[2] - atoms[metal_indices[-3]].position[2])
            oxygen_positions = [atoms[metal_indices[-2]].position + bond_vector,
                                atoms[metal_indices[-1]].position + bond_vector]
            hydrogen_positions = [pos + (1.0, 0.0, 0.0) for pos in oxygen_positions]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o.json', 'restart-oh.json')

        elif '7_Pyramidal_LT' in dir:
            bond_vector = atoms[oxygen_indices[-1]].position - atoms[metal_indices[-1]].position
            oxygen_positions = [atoms[metal_indices[-1]].position + (0.0, 0.0, np.linalg.norm(bond_vector))]
            hydrogen_positions = [oxygen_positions[0] + (0.0, 1.0, 0.0)]
            add_atoms_and_save(atoms, oxygen_positions, hydrogen_positions, bond_vector[2], 'restart-o.json', 'restart-oh.json')

        print(f'Saved json files in {dir}')
    except Exception as e:
        print(f"Error processing directory {dir}: {e}")
