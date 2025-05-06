import os
import glob
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.geometry import get_duplicate_atoms
from ase.build import make_supercell, surface
from ase.build.tools import sort
import socket

# 서버 주소 가져오기
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/'
elif user_name == 'hailey':
    root = '/Users/hailey/Desktop/7_V_bulk'
else:
    raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

coords = [
    {'coord_dir': '1_Tetrahedral_WZ',  'P': np.array([[0,0,1],[1,2,0],[-1,0,0]]), 'n_layers': 5, 'm_fix': 0},
    {'coord_dir': '2_Tetrahedral_ZB',  'P': np.array([[0.5,-0.5,-0.5],[1.5,1.5,-2.5],[1.5,-1.5,0.5]])*2, 'n_layers': 2, 'm_fix': 0},
]

rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

for coord in coords:
    coord_dir = coord['coord_dir']
    P = coord['P']
    facet = (0,0,1)
    n_layers = coord['n_layers']
    fixed_atom_count = coord['m_fix']
    vacuum = 9.0 # Å
    for row in rows:
        metals = rows[row]
        for m, metal in enumerate(metals):
            numb = str(m).zfill(2)
            path = os.path.join(root, coord_dir, row, f'{numb}_{metal}')
            # unmatched_path = os.path.join(path, 'unmatched')
            json_path = os.path.join(path, 'final_with_calculator.json')
            if os.path.exists(json_path): # and not os.path.exists(unmatched_path) :
                # Step 1: Read the original structure
                atoms = read(json_path)

                # Step 2: Apply the supercell transformation
                atoms = make_supercell(atoms, P)

                # Step 3: Create a slab with the (010) facet
                slab = surface(atoms, indices=facet, layers=n_layers, vacuum=vacuum, periodic=True)

                # Step 4: Center the slab along the z-direction
                slab.center(axis=2)
                slab.wrap()
                slab = sort(slab)

                # get_duplicate_atoms(slab, cutoff=0.1, delete=True)
                # sorted_atoms = sorted(slab, key=lambda atom: atom.position[2])
                # indices_to_fix = [atom.index for atom in sorted_atoms[:fixed_atom_count]]
                # slab.set_constraint(FixAtoms(indices=indices_to_fix))

                # Step 5: Save the final slab structure
                output_path = os.path.join(path, "slab.traj")
                write(output_path, slab)

                print(f"Created slab.traj at: {coord_dir}/{numb}_{metal}")