#!/bin/env python

import os
from ase.io import read, write

root = '/pscratch/sd/j/jiuy97/7_V_bulk'
coords = {
    'WZ': '1_Tetrahedral_WZ', 
    'ZB': '2_Tetrahedral_ZB', 
    'TN': '3_SquarePlanar_TN', 
    'PD': '4_SquarePlanar_PD', 
    'NB': '5_SquarePlanar_NB', 
    'RS': '6_Octahedral_RS',
    'LT': '7_Pyramidal_LT', 
    'AQ': '8_Tetrahedral_AQ', 
    'AU': '9_SquarePlanar_AU'
}
metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

def main():
    for coord in coords.keys():
        for row in metals.keys():
            for m, metal in enumerate(metals[row]):
                numb = str(m).zfill(2)
                dir_path = os.path.join(root, coords[coord], row, numb+'_'+metal)
                atoms_path = os.path.join(dir_path, 'final_with_calculator.json')
                if os.path.exist(atoms_path):
                    atoms = read(atoms_path)
                print(dir_path)

if __name__ == "__main__":
    main()
