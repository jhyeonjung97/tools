from ase.io import read, write
import numpy as np
import os

# root = '/Users/hailey/Desktop/1_cation/3_PES/1_cation_water/0_'
root = '/sdf/data/suncat/suncatlustre/hailey/1_cation/3_PES/1_cation_water/0_'
dirs = ['1_Li', '2_Na', '3_K', '4_Rb', '5_Cs']


for i, dir in enumerate(dirs):

    cation_file = os.path.join(root, dir, 'restart.json')
    cation_atoms = read(cation_file)

    target_index = 19
    z_displacement = cation_atoms.cell.lengths()[2] / 6 * 2
    for atom in cation_atoms:
        if atom.symbol == 'Li' or atom.symbol == 'Na':
            target_index = 18
            z_displacement = cation_atoms.cell.lengths()[2] / 6 * 3        

    cation_atoms.positions[:, 2] -= z_displacement
    cation_atoms.wrap()

    atoms = read('restart.json')
    
    l1 = atoms.cell.lengths()[0]
    l2 = atoms.cell.lengths()[1]
    l3 = atoms.cell.lengths()[2]
    a1 = atoms.cell.angles()[0]
    a2 = atoms.cell.angles()[1]
    a3 = atoms.cell.angles()[2]
    atoms.cell = (l1, l2, 40, a1, a2, a3)

    center_x = l1/4
    center_y = l2/4*np.sqrt(3)
    shift_x = center_x - atoms.positions[target_index, 0]
    shift_y = center_y - atoms.positions[target_index, 1]
    atoms.positions[:, 0] += shift_x
    atoms.positions[:, 1] += shift_y

    positions = atoms.get_positions()
    z_sorted_indices = np.argsort(positions[:, 2])
    indices_to_delete = z_sorted_indices[:9]

    del atoms[indices_to_delete]
    atoms.wrap()

    min_z = atoms.positions[:, 2].min()
    shift = 0.1 - min_z
    atoms.positions[:, 2] += shift
    max_z = atoms.positions[:, 2].max()

    l3 = cation_atoms.cell[2][2]
    cation_atoms.set_cell([l1, l2, l3, a1, a2, a3], scale_atoms=True)

    shift = max_z + 1.2 - cation_atoms.positions[:, 2].min()
    cation_atoms.positions[:, 2] += shift
    new_atoms = atoms + cation_atoms

    write(f'restart_hol{i+1}.json', new_atoms)