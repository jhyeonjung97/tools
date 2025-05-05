import os
import glob
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.geometry import get_duplicate_atoms
from ase.build import make_supercell, surface
from ase.build.tools import sort
import socket

# Define root directory and coordination type
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/7_V_bulk'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")

coord = '6_Octahedral_RS'

rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

for row in rows:
    metals = rows[row]
    for m, metal in enumerate(metals):
        numb = str(m).zfill(2)
        path = os.path.join(root, coord, row, f'{numb}_{metal}')
        atoms_path = os.path.join(root, coord, row, f'a{numb}.vasp')
        
        if not os.path.exists(atoms_path):
            continue
            
        atoms = read(atoms_path)
        
        l1 = atoms.cell.lengths()[0]
        l2 = atoms.cell.lengths()[1]
        l3 = atoms.cell.lengths()[2]
        a1 = atoms.cell.angles()[0]
        a2 = atoms.cell.angles()[1]
        a3 = atoms.cell.angles()[2]
        atoms.cell = (l1/2, l2/2, l3, a1, a2, a3)
        
        atoms.wrap()
        get_duplicate_atoms(atoms, cutoff=1.0, delete=True)
        
        atoms.positions[:, 0] += l1/3
        atoms.positions[:, 1] += l2/3
        atoms.wrap()
        get_duplicate_atoms(atoms, cutoff=1.0, delete=True)
        
        z_positions = atoms.positions[:, 2]
        min_z = z_positions.min()
        max_z = z_positions.max()
        slab_thickness = max_z - min_z
        
        new_height = slab_thickness + 18
        cell = atoms.cell.array.copy()
        cell[2, :] = [0.0, 0.0, new_height]
        atoms.set_cell(cell)
        atoms.center(axis=2)
        atoms = sort(atoms)
        
        sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
        indices_to_fix = [atom.index for atom in sorted_atoms[:8]]
        atoms.set_constraint(FixAtoms(indices=indices_to_fix))
        
        if sorted_atoms[0].symbol == 'O':
            atoms.positions[:, 0] -= sorted_atoms[0].position[0]
            atoms.positions[:, 1] -= sorted_atoms[0].position[1]
        elif sorted_atoms[1].symbol == 'O':
            atoms.positions[:, 0] -= sorted_atoms[1].position[0]
            atoms.positions[:, 1] -= sorted_atoms[1].position[1]
        elif sorted_atoms[2].symbol == 'O':
            atoms.positions[:, 0] -= sorted_atoms[2].position[0]
            atoms.positions[:, 1] -= sorted_atoms[2].position[1]
        atoms.wrap()
        
        output_path = os.path.join(path, "slab.traj")
        write(output_path, atoms)

        # print(f"Created slab.traj at: 2_Tetrahedral_ZB/{numb}_{metal}")
        
        if len(atoms) != 24:
            print(f"Created slab.traj at: {coord}/{row}/{numb}_{metal}")
