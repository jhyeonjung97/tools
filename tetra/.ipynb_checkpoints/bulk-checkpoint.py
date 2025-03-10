#!/bin/env python

import os
import numpy as np
import pandas as pd
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

df = pd.DataFrame() # Madelung, GP, hexa
def main():
    for coord in coords.keys():
        CN = 4 if coord != 'RS' else 6
        if coord == 'AU':
            ON = 1.5
        elif coord == 'AQ':
            ON = 2
        else:
            ON = 1
        for row in metals.keys():
            for m, metal in enumerate(metals[row]):
                numb = str(m).zfill(2)
                item = coord+row+metal
                df.loc[item, ['coord', 'CN', 'ON', 'row', 'metal']] = coord, CN, ON, row, metal
                
                dir_path = os.path.join(root, coords[coord], row, numb+'_'+metal)
            
                atoms_path = os.path.join(dir_path, 'isif2/final_with_calculator.json')
                if os.path.exists(atoms_path):
                    atoms = read(atoms_path)
                    energy = atoms.get_total_energy()
                    volume = atoms.get_volume()
                    df.loc[item, ['energy', 'volume']] = energy, volume
            
                chg_path = os.path.join(dir_path, 'isif2/atoms_bader_charge.json')
                if os.path.exists(chg_path):
                    atoms = read(chg_path)
                    chgs = atoms.get_charges()
                    chg = np.mean([chgs[atom.index] for atom in atoms if atom.symbol == metal])
                    df.loc[item, 'chg'] = chg
            
                mag_path = os.path.join(dir_path, 'isif2/moments.json')
                if os.path.exists(mag_path):
                    atoms = read(mag_path)
                    mags = atoms.get_magnetic_moments()
                    mag = np.mean([mags[atom.index] for atom in atoms if atom.symbol == metal])
                    df.loc[item, 'mag'] = mag
            
                icohp_path = os.path.join(dir_path, 'icohp.txt')
                icobi_path = os.path.join(dir_path, 'icobi.txt')
                icoop_path = os.path.join(dir_path, 'icoop.txt')
                if os.path.exists(icohp_path):
                    icohp, bond = parse_icohp(icohp_path)
                    icobi, bond = parse_icohp(icobi_path)
                    icoop, bond = parse_icohp(icoop_path)
                    print(icohp, bond)


def parse_icohp(file_path):
    distances = []
    icohps = []

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7 or parts[0] == "label":  
                # Skip header and invalid lines
                continue  

            try:
                icohp = float(parts[-2])  # -ICOHP is the second last value
                distance = float(parts[-1])  # Distance is the last value
                icohps.append(icohp)
                distances.append(distance)
            except ValueError:
                continue  # Skip lines that can't be converted to float

    if icohps and distances:
        avg_icohp = np.mean(icohps)
        avg_distance = np.mean(distances)
    else:
        avg_icohp = avg_distance = None  # Return None if no valid data

    return avg_icohp, avg_distance

            
if __name__ == "__main__":
    main()
