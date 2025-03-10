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

# df = pd.DataFrame() # Madelung, GP, hexa
# df = pd.DataFrame(columns=['coord', 'CN', 'ON', 'row', 'metal', 'energy', 'volume', 'chg', 'mag', 'icohp', 'icobi', 'icoop'])
# df = df.astype({'coord': 'string', 'row': 'string', 'metal': 'string', 'CN': 'float64', 'ON': 'float64'})
df = pd.DataFrame(columns=['coord', 'CN', 'ON', 'row', 'metal', 
                           'energy', 'volume', 'chg', 'mag', 
                           'bond', 'icohp', 'icobi', 'icoop', 'madelung', 'grosspop'],
                  dtype='object')

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
                    chgs = atoms.get_initial_charges()
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
                madelung_path = os.path.join(dir_path, 'MadelungEnergies.lobster')
                grosspop_path = os.path.join(dir_path, 'GROSSPOP.lobster')
                if os.path.exists(icohp_path) and os.path.getsize(icohp_path) != 0:
                    icohp, bond = parse_icohp(icohp_path)
                    icobi = parse_icohp(icobi_path)
                    icoop = parse_icohp(icoop_path)
                    madelung = parse_madelung(madelung_path)
                    grosspop = parse_grosspop(grosspop_path, metal)
                    df.loc[item, ['bond', 'icohp', 'icobi', 'icoop']] = bond, icohp, icobi, icoop
                    df.loc[item, ['madelung', 'grosspop']] = madelung, grosspop
                    print(dir_path)
                    print(madelung, grosspop)
    
def parse_icohp(file_path):
    distances, icohps = [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7 or parts[0] == "label":
                continue  
            icohps.append(float(parts[-2]))  
            distances.append(float(parts[-1]))
    avg_icohp = np.mean(icohps) if icohps else None
    avg_distance = np.mean(distances) if distances else None
    return avg_icohp, avg_distance

def parse_madelung(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3 and parts[0].replace('.', '', 1).isdigit():
                madelung = float(parts[2])  # Third column is Loewdin energy
                return madelung
    return None
                                      
import numpy as np

def parse_grosspop(file_path, metal):
    elements, loewdin_totals = [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 5 and (parts[1] == metal or parts[1] == 'O'):
                elements.append(parts[1])  
            if len(parts) == 5 and parts[2] == 'total':
                loewdin_totals.append(float(parts[-1]))
    print(elements)
    print(loewdin_totals)
    loewdin_gp = [loewdin_totals[i] for i in range(len(elements)) if elements[i] == metal]
    return np.mean(loewdin_gp) if loewdin_gp else None

    
if __name__ == "__main__":
    main()
