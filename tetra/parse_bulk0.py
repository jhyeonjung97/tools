#!/bin/env python

import os
import socket
import numpy as np
import numpy.ma as ma
import pandas as pd
from ase.io import read, write
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mendeleev import element

# 서버 주소 가져오기
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk/'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/7_V_bulk/'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")
save_path = os.path.join(root, 'figures')

coords_data = [
    {'coord': 'WZ', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '1_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
    {'coord': 'ZB', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '2_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
    {'coord': 'TN', 'CN': 4, 'ON': 2, 'MN': 4, 'coord_dir': '3_SquarePlanar_TN', 'zorder': 3, 'marker': 'o', 'color': 'dodgerblue',},
    {'coord': 'PD', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '4_SquarePlanar_PD', 'zorder': 2, 'marker': 'o', 'color': 'deepskyblue',},
    {'coord': 'NB', 'CN': 4, 'ON': 2, 'MN': 6, 'coord_dir': '5_SquarePlanar_NB', 'zorder': 1, 'marker': 's', 'color': 'limegreen',},
    {'coord': 'RS', 'CN': 6, 'ON': 2, 'MN': 2, 'coord_dir': '6_Octahedral_RS',   'zorder': 6, 'marker': 'd', 'color': 'orchid',},
    {'coord': 'LT', 'CN': 4, 'ON': 2, 'MN': 2, 'coord_dir': '7_Pyramidal_LT',    'zorder': 0, 'marker': 'h', 'color': 'silver',},
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
    {'column': 'coord', 'png_name': 'coordination', 'ylabel': 'Coordination'},
    {'column': 'row',   'png_name': 'row',          'ylabel': 'Row'},
    {'column': 'numb',  'png_name': 'number',       'ylabel': 'Number'},
    {'column': 'metal', 'png_name': 'metal',        'ylabel': 'Metal'},
    {'column': 'CN',    'png_name': 'coordination_number', 'ylabel': 'Coordination Number'},
    {'column': 'ON',    'png_name': 'oxidation_number',    'ylabel': 'Oxidation Number'},
    {'column': 'cell',  'png_name': 'cell',         'ylabel': 'Cell Distortion'},
]
columns = pd.DataFrame(columns_data).set_index('column')
columns.index.name = None

df = pd.DataFrame(columns=columns.index, dtype='object')
int_cols = ['CN', 'ON']
str_cols = ['coord', 'row', 'numb', 'metal']
float_cols = ['cell']

def main():
    global df
    
    for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']:
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
                    
                    if coord in ['WZ', 'TN', 'PD', 'LT']:
                        a = atoms.cell.cellpar()[0]
                        c = atoms.cell.cellpar()[2]
                        df.loc[item, 'cell'] = c/a
                    elif coord in ['ZB', 'RS']:
                        df.loc[item, 'cell'] = atoms.cell.cellpar()[3]
                    elif coord in ['NB']:
                        df.loc[item, 'cell'] = atoms.cell.cellpar()[3]
                    
                df.to_csv(f'{save_path}/bulk_data.csv', sep=',')
                df[float_cols] = df[float_cols].astype(float).round(2)
                df.to_csv(f'{save_path}/bulk_data.tsv', sep='\t', float_format='%.2f')
    
    print(df)
    plot_by_metal_row(df, save_path)
    plot_by_coordination(df, save_path)
    
def plot_by_metal_row(df, save_path):
    for row in ['fm', '3d', '4d', '5d']:
        plt.figure(figsize=(4, 3))
        for c, coord in enumerate(['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']):
            zorder=coords.loc[coord, 'zorder']
            marker=coords.loc[coord, 'marker']
            color=coords.loc[coord, 'color']
            subset = df[(df['coord'] == coord) & (df['row'] == row)]
            plt.plot(subset['numb'], subset['cell'], marker=marker, color=color, 
                     linestyle='-', label=coord, zorder=zorder)
                            
        plt.xticks(np.arange(len(metals[row])), metals[row])
        plt.xlabel("Metal Index")
        plt.ylabel(columns.loc['cell', 'ylabel'])
        plt.legend()
        plt.tight_layout()
        png_name = f"bulk_{row}_{columns.loc['cell', 'png_name']}.png"
        plt.savefig(f"{save_path}/{png_name}")
        plt.close()
        print(f"Figure saved as {png_name}")
            
def plot_by_coordination(df, save_path):        
    for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']:            
        marker=coords.loc[coord, 'marker']
        base_color = coords.loc[coord, 'color']
        cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
        colors = cmap(np.linspace(0.0, 0.6, 3))
        plt.figure(figsize=(4, 3))
        for r, row in enumerate(['3d', '4d', '5d']):
            color = 'lightgray' if row == 'fm' else colors[r]
            subset = df[(df['coord'] == coord) & (df['row'] == row)]           
            plt.plot(subset['numb'], subset['cell'], marker=marker, color=color, 
                     linestyle='-', label=row)
                    
        plt.xticks(np.arange(len(indice)), indice)
        plt.xlabel("Metal Index")
        plt.ylabel(columns.loc['cell', 'ylabel'])
        plt.legend()
        plt.tight_layout()
        png_name = f"bulk_{coord}_{columns.loc['cell', 'png_name']}.png"
        plt.savefig(f"{save_path}/{png_name}")
        plt.close()
        print(f"Figure saved as {png_name}")

if __name__ == "__main__":
    main()