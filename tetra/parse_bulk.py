#!/bin/env python

import os
import numpy as np
import pandas as pd
from ase.io import read, write
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

root = '/pscratch/sd/j/jiuy97/7_V_bulk'
save_path = os.path.join(root, 'figures')

coords_data = [
    {'coord': 'WZ', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '1_Tetrahedral_WZ',  'marker': '>', 'color': 'tab:red',   },
    {'coord': 'ZB', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '2_Tetrahedral_ZB',  'marker': '<', 'color': 'tab:orange',},
    {'coord': 'TN', 'CN': 4, 'ON': 2, 'MN': 4, 'coord_dir': '3_SquarePlanar_TN', 'marker': 'o', 'color': 'tab:olive', },
    {'coord': 'PD', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '4_SquarePlanar_PD', 'marker': 's', 'color': 'tab:green', },
    {'coord': 'NB', 'CN': 4, 'ON': 2, 'MN': 6, 'coord_dir': '5_SquarePlanar_NB', 'marker': 'p', 'color': 'tab:blue',  },
    {'coord': 'RS', 'CN': 6, 'ON': 2, 'MN': 2, 'coord_dir': '6_Octahedral_RS',   'marker': 'd', 'color': 'tab:purple',},
    {'coord': 'LT', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '7_Pyramidal_LT',    'marker': 'h', 'color': 'tab:brown', },
    {'coord': 'AQ', 'CN': 4, 'ON': 4, 'MN': 6, 'coord_dir': '8_Tetrahedral_AQ',  'marker': '^', 'color': 'tab:pink',  },
    {'coord': 'AU', 'CN': 4, 'ON': 3, 'MN': 4, 'coord_dir': '9_SquarePlanar_AU', 'marker': 'v', 'color': 'tab:cyan',  },
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

columns_data = [
    {'column': 'coord',    'png_name': 'coordination',        'ylabel': 'Coordination'},
    {'column': 'row',      'png_name': 'row',                 'ylabel': 'Row'},
    {'column': 'numb',     'png_name': 'number',              'ylabel': 'Number'},
    {'column': 'metal',    'png_name': 'metal',               'ylabel': 'Metal'},
    {'column': 'CN',       'png_name': 'coordination_number', 'ylabel': 'Coordination Number'},
    {'column': 'ON',       'png_name': 'oxidation_number',    'ylabel': 'Oxidation Number'},
    {'column': 'energy',   'png_name': 'energy',              'ylabel': 'Energy (eV)'},
    {'column': 'form',     'png_name': 'formation_energy',    'ylabel': 'Formation Energy (eV)'},
    {'column': 'volume',   'png_name': 'volume',              'ylabel': 'Volume (Å³)'},
    {'column': 'cell',     'png_name': 'cell',                'ylabel': 'Cell'},
    {'column': 'chg',      'png_name': 'bader_charge',        'ylabel': 'Bader Charge (e⁻)'},
    {'column': 'mag',      'png_name': 'magnetic_moments',    'ylabel': 'Magnetic Moments (μB)'},
    {'column': 'l_bond',   'png_name': 'bond_length',         'ylabel': 'Bond Length (Å)'},
    {'column': 'n_bond',   'png_name': 'number_of_bonds',     'ylabel': 'Number of Bonds per Metal'},
    {'column': '-ICOHPm',  'png_name': 'icohp_per_metal',     'ylabel': '-ICOHP per Metal (eV)'},
    {'column': 'ICOBIm',   'png_name': 'icobi_per_metal',     'ylabel': 'ICOBI per Metal'},
    {'column': '-ICOOPm',  'png_name': 'icoop_per_metal',     'ylabel': '-ICOOP per Metal (eV)'},
    {'column': '-ICOHPn',  'png_name': 'icohp_per_bond',      'ylabel': '-ICOHP per Bond (eV)'},
    {'column': 'ICOBIn',   'png_name': 'icobi_per_bond',      'ylabel': 'ICOBI per Bond'},
    {'column': '-ICOOPn',  'png_name': 'icoop_per_bond',      'ylabel': '-ICOOP per Bond (eV)'},
    {'column': 'madelung', 'png_name': 'madelung',            'ylabel': 'Madelung Energy (Loewdin, eV)'},
    {'column': 'grosspop', 'png_name': 'gross_population',    'ylabel': 'Gross Population (Loewdin, e⁻)'},
]
columns = pd.DataFrame(columns_data).set_index('column')
columns.index.name = None

df = pd.DataFrame(columns=columns.index, dtype='object')
int_cols = ['CN', 'ON', 'n_bond']
float_cols = ['energy', 'form', 'volume', 'cell', 'chg', 'mag', 'l_bond', 'n_bond', '-ICOHPm', 'ICOBIm', '-ICOOPm', '-ICOHPn', 'ICOBIn', '-ICOOPn', 'madelung', 'grosspop']

metal_df = pd.read_csv('~/bin/tools/tetra/metal-data.tsv', sep='\t', index_col=0)

E_H2O = -14.23919983
E_H2 = -6.77409008

ZPE_H2O = 0.558
ZPE_H2 = 0.273

Cp_H2O = 0.10
Cp_H2 = 0.09

Ref_H2 = E_H2 + ZPE_H2 + Cp_H2
Ref_H2O = E_H2O + ZPE_H2O + Cp_H2O
Ref_O2 = 2*Ref_H2O - 2*Ref_H2 + 4.92
Ref_O = Ref_H2O - Ref_H2

def main():
    global df
    
    if os.path.exists(f'{save_path}/bulk_data.csv'):
        df = pd.read_csv(f'{save_path}/bulk_data.csv')
    else:
        # for coord in coords.index:
        for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS']:
            CN = coords.loc[coord, 'CN']
            ON = coords.loc[coord, 'ON']
            MN = coords.loc[coord, 'MN']
            coord_dir = coords.loc[coord, 'coord_dir']
            
            for row in metals.keys():
                for m, metal in enumerate(metals[row]):
                    numb = str(m).zfill(2)
                    item = coord+row+numb
                    df.loc[item, ['coord', 'row', 'numb', 'metal', 'CN', 'ON']] = coord, row, m, metal, CN, ON 
                    dir_path = os.path.join(root, coord_dir, row, numb+'_'+metal)
                    
                    atoms_path = os.path.join(dir_path, 'isif2/final_with_calculator.json')                
                    if os.path.exists(atoms_path):
                        atoms = read(atoms_path)
                        energy = atoms.get_total_energy()
                        volume = atoms.get_volume()
                        df.loc[item, ['energy', 'volume']] = energy/MN, volume/MN

                        formation = energy/MN - metal_df.loc[metal, 'E'] - (Ref_O2 / 2) * (ON /2)
                        df.loc[item, 'form'] = formation

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
    plot_by_metal_row(df, save_path)
    plot_by_coordination(df, save_path)
    
def plot_by_metal_row(df, save_path):
    for row in metals.keys():
        for col in columns.index:
            plt.figure(figsize=(8, 6))
            for coord in coords.index:
                marker=coords.loc[coord, 'marker']
                color=coords.loc[coord, 'color']
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, linestyle='-', label=coord)
                for m, metal in enumerate(metals[row]):
                    numb = str(m).zfill(2)
                    print(numb, metal_df.loc[metal, 'Eform'])
                    # plt.scatter(m, metal_df.loc[metal, 'Eform'], marker=marker,
                    #             edgecolors=color, facecolors='white', label='exp.')
                        
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(columns.loc[col, 'ylabel'])
            plt.legend()
            plt.tight_layout()
            png_name = f"bulk_{row}_{columns.loc[col, 'png_name']}.png"
            plt.savefig(f"{save_path}/{png_name}")
            plt.close()
            print(f"Figure saved as {png_name}")
            
def plot_by_coordination(df, save_path):        
    for coord in coords.index:
        for col in columns.index:
            marker=coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', ['white', base_color])
            colors = cmap(np.linspace(0.3, 1, 3))
            plt.figure(figsize=(8, 6))
            for r, row in enumerate(['fm', '3d', '4d', '5d']):
                color = 'lightgray' if row == 'fm' else colors[r-1]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, linestyle='-', label=row)
                for m, metal in enumerate(metals[row]):
                    numb = str(m).zfill(2)
                    print(numb, metal_df.loc[metal, 'Eform'])
                    # plt.scatter(m, metal_df.loc[metal, 'Eform'], marker=marker,
                    #             edgecolors=color, facecolors='white', label='exp.')
                        
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(columns.loc[col, 'ylabel'])
            plt.legend()
            plt.tight_layout()
            png_name = f"bulk_{coord}_{columns.loc[col, 'png_name']}.png"
            plt.savefig(f"{save_path}/{png_name}")
            plt.close()
            print(f"Figure saved as {png_name}")

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
