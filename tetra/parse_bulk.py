#!/bin/env python

import os
import numpy as np
import pandas as pd
from ase.io import read, write
import matplotlib.pyplot as plt

root = '/pscratch/sd/j/jiuy97/7_V_bulk'
save_path = os.path.join(root, 'figures')

coords_data = [
    {'coord': 'WZ', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '1_Tetrahedral_WZ'},
    {'coord': 'ZB', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '2_Tetrahedral_ZB'},
    {'coord': 'TN', 'CN': 4, 'ON': 2, 'MN': 4, 'coord_dir': '3_SquarePlanar_TN'},
    {'coord': 'PD', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '4_SquarePlanar_PD'},
    {'coord': 'NB', 'CN': 4, 'ON': 2, 'MN': 6, 'coord_dir': '5_SquarePlanar_NB'},
    {'coord': 'RS', 'CN': 6, 'ON': 2, 'MN': 2, 'coord_dir': '6_Octahedral_RS'},
    {'coord': 'LT', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '7_Pyramidal_LT'},
    {'coord': 'AQ', 'CN': 4, 'ON': 4, 'MN': 6, 'coord_dir': '8_Tetrahedral_AQ'},
    {'coord': 'AU', 'CN': 4, 'ON': 3, 'MN': 4, 'coord_dir': '9_SquarePlanar_AU'}
]
coords = pd.DataFrame(coords_data).set_index('coord')
coords.index.name = None

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])]

columns={
    'coord': 'Coordination', 
     'row': 'Row', 
     'numb': 'Number', 
     'metal': 'Metal', 
     'CN': 'Coordination Number', 
     'ON': 'Oxidation Number', 
     'energy': 'Energy (eV)', 
     'volume': 'Volume (Å³)', 
     'cell': 'Cell', 
     'chg': 'Bader Charge (e⁻)', 
     'mag': 'Magnetic Moments (μB)',
     'l_bond': 'Bond Length (Å)', 
     'n_bond': 'Number of Bonds per Metal', 
     '-ICOHPm': '-ICOHP per Metal (eV)', 
     'ICOBIm': 'ICOBI per Metal', 
     '-ICOOPm': '-ICOOP per Metal (eV)', 
     '-ICOHPn': '-ICOHP per Bond (eV)', 
     'ICOBIn': 'ICOBI per Metal', 
     '-ICOOPn': '-ICOOP per Bond (eV)', 
     'madelung': 'Madelung Energy (Loewdin, eV)',
     'grosspop': 'Gross Population (Loewdin, e⁻)',
}

df = pd.DataFrame(columns=['coord', 'row', 'numb', 'metal', 'CN', 'ON', 
                           'energy', 'volume', 'cell', 'chg', 'mag',
                           'l_bond', 'n_bond', '-ICOHPm', 'ICOBIm', '-ICOOPm', '-ICOHPn', 'ICOBIn', '-ICOOPn', 'madelung', 'grosspop'],
                  dtype='object')
int_cols = ['CN', 'ON', 'n_bond']
float_cols = ['energy', 'volume', 'cell', 'chg', 'mag', 'l_bond', 'n_bond', '-ICOHPm', 'ICOBIm', '-ICOOPm', '-ICOHPn', 'ICOBIn', '-ICOOPn', 'madelung', 'grosspop']

def main():
    global df, coords
    
    for coord in coords.index:
        CN = coords.loc[coord, 'CN']
        ON = coords.loc[coord, 'ON']
        MN = coords.loc[coord, 'MN']
        coord_dir = coords.loc[coord, 'coord_dir']

        if os.path.exists(f'{save_path}/bulk_data.csv'):
            df = pd.read_csv(f'{save_path}/bulk_data.csv')
            break

        for row in metals.keys():
            for m, metal in enumerate(metals[row]):
                numb = str(m).zfill(2)
                item = coord+row+numb
                df.loc[item, ['coord', 'row', 'numb', 'metal', 'CN', 'ON']] = coord, row, numb, metal, CN, ON 
                dir_path = os.path.join(root, coord_dir, row, numb+'_'+metal)
                
                atoms_path = os.path.join(dir_path, 'isif2/final_with_calculator.json')                
                if os.path.exists(atoms_path):
                    atoms = read(atoms_path)
                    energy = atoms.get_total_energy()
                    volume = atoms.get_volume()
                    df.loc[item, ['energy', 'volume']] = energy/MN, volume/MN

                    if coord in ['WZ', 'TN', 'PD', 'LT', 'AQ']:
                        a = atoms.cell.cellpar()[0]
                        c = atoms.cell.cellpar()[2]
                        df.loc[item, 'cell'] = c/a
                    elif coord in ['ZB', 'NB', 'RS']:
                        df.loc[item, 'cell'] = atoms.cell.cellpar()[3]
                
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
                    mag = np.mean([abs(mags[atom.index]) for atom in atoms if atom.symbol == metal])
                    df.loc[item, 'mag'] = mag
            
                icohp_path = os.path.join(dir_path, 'icohp.txt')
                icobi_path = os.path.join(dir_path, 'icobi.txt')
                icoop_path = os.path.join(dir_path, 'icoop.txt')
                madelung_path = os.path.join(dir_path, 'MadelungEnergies.lobster')
                grosspop_path = os.path.join(dir_path, 'GROSSPOP.lobster')
                if os.path.exists(icohp_path) and os.path.getsize(icohp_path) != 0:
                    icohp, bond, nbond = parse_icohp(icohp_path)
                    icobi, _, _ = parse_icohp(icobi_path)
                    icoop, _, _ = parse_icohp(icoop_path)
                    madelung = parse_madelung(madelung_path)
                    grosspop = parse_grosspop(grosspop_path, metal)
                    df.loc[item, ['l_bond', 'n_bond', '-ICOHPn', 'ICOBIn', '-ICOOPn', 'madelung', 'grosspop']] = bond, nbond, icohp, icobi, icoop, madelung/MN, grosspop
                    df.loc[item, ['-ICOHPm', 'ICOBIm', '-ICOOPm']] = icohp*nbond, icobi*nbond, icoop*nbond
                    # if CN != nbond:
                    #     print(dir_path)

                df.to_csv(f'{save_path}/bulk_data.csv', sep=',')
                # df[int_cols] = df[int_cols].astype(int)
                df[float_cols] = df[float_cols].astype(float).round(2)
                df.to_csv(f'{save_path}/bulk_data.tsv', sep='\t', float_format='%.2f')
                
    print(df)
    plot_by_metal_row(df, coords, save_path)
    plot_by_coordination(df, coords, save_path)
        
def plot_by_metal_row(df, save_path):        
    for coord in coords.index:
        for col in df.columns:
            plt.figure(figsize=(8, 6))
            for row in metals.keys():
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker='o', linestyle='-', label=row)
                
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(col.title())
            plt.legend()
            plt.tight_layout()
            # plt.savefig(f"{save_path}/bulk_{coord}_{col}.png")
            plt.savefig(f"{save_path}/bulk_{coord}.png")
            plt.close()
            
def plot_by_coordination(df, save_path):
    for row in metals.keys():
        for col in df.columns:
            plt.figure(figsize=(8, 6))
                for coord in coords.index:
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker='o', linestyle='-', label=row)
                
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(col.title())
            plt.legend()
            plt.tight_layout()
            # plt.savefig(f"{save_path}/bulk_{row}_{col}.png")
            plt.savefig(f"{save_path}/bulk_{row}.png")
            plt.close()
            
def parse_icohp(file_path):
    distances, icohps = [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 8 and parts[0] != 'label':
                icohps.append(float(parts[-2]))  
                distances.append(float(parts[-1]))
    nbond = len(icohps) if icohps else np.nan
    avg_icohp = np.mean(icohps) if icohps else np.nan
    avg_distance = np.mean(distances) if distances else np.nan
    return avg_icohp, avg_distance, nbond

def parse_madelung(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3 and parts[0].replace('.', '', 1).isdigit():
                madelung = float(parts[2])
                return madelung
    return np.nan
                                      
def parse_grosspop(file_path, metal):
    elements, loewdin_totals = [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 5 and (parts[1] == metal or parts[1] == 'O'):
                elements.append(parts[1])  
            if len(parts) == 3 and parts[0] == 'total':
                loewdin_totals.append(float(parts[-1]))
    loewdin_gp = [loewdin_totals[i] for i in range(len(elements)) if elements[i] == metal]
    return np.mean(loewdin_gp) if loewdin_gp else np.nan

    
if __name__ == "__main__":
    main()
