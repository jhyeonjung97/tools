#!/bin/env python

import os
import socket
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mendeleev import element

hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk/figures'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/7_V_bulk/figures'
else:
    raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

data_path = f'{root}/bulk_data_cfse.csv'
save_path = f'{root}'

# 데이터 읽기
df = pd.read_csv(data_path, index_col=0)

# 좌표 정보 설정
coords_data = [
    {'coord': 'WZ', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '1_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'orangered'},
    {'coord': 'ZB', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '2_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'darkorange'},
    {'coord': 'TN', 'CN': 4, 'OS': 2, 'MN': 4, 'coord_dir': '3_SquarePlanar_TN', 'zorder': 3, 'marker': 's', 'color': 'dodgerblue'},
    {'coord': 'PD', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '4_SquarePlanar_PD', 'zorder': 2, 'marker': 's', 'color': 'deepskyblue'},
    {'coord': 'NB', 'CN': 4, 'OS': 2, 'MN': 6, 'coord_dir': '5_SquarePlanar_NB', 'zorder': 1, 'marker': 's', 'color': 'limegreen'},
    {'coord': 'RS', 'CN': 6, 'OS': 2, 'MN': 2, 'coord_dir': '6_Octahedral_RS',   'zorder': 6, 'marker': 'd', 'color': 'orchid'},
    {'coord': 'LT', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '7_Pyramidal_LT',    'zorder': 0, 'marker': 'o', 'color': 'gold'},
    {'coord': '+3', 'CN': 6, 'OS': 3, 'MN': 4, 'coord_dir': '1_Octahedral_+3',   'zorder': 1, 'marker': 'd', 'color': (0.9, 0.9, 0.9)},
    {'coord': '+4', 'CN': 6, 'OS': 4, 'MN': 2, 'coord_dir': '2_Octahedral_+4',   'zorder': 2, 'marker': 'd', 'color': (0.8, 0.8, 0.8)},
    {'coord': '+5', 'CN': 6, 'OS': 5, 'MN': 4, 'coord_dir': '3_Octahedral_+5',   'zorder': 3, 'marker': 'd', 'color': (0.7, 0.7, 0.7)},
    {'coord': '+6', 'CN': 6, 'OS': 6, 'MN': 1, 'coord_dir': '4_Octahedral_+6',   'zorder': 4, 'marker': 'd', 'color': (0.6, 0.6, 0.6)},
]

coords = pd.DataFrame(coords_data).set_index('coord')
coords.index.name = None

# 금속 정보 설정
metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
}
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])]

# 에너지 계산 함수
def calculate_energies(df):
    # 새로운 컬럼 추가
    df['sublimation_energy'] = 0.0  # kJ/mol -> eV
    df['ionization_energy'] = 0.0   # kJ/mol -> eV
    df['lattice_energy'] = 0.0      # eV
    
    # 산소 원자의 원자화 에너지 (eV)
    oxygen_atomization = 2.58  # eV
    
    # 산소의 전자 친화도 (eV)
    oxygen_electron_affinity = 1.46  # eV
    
    for idx, row in df.iterrows():
        metal = row['metal']
        oxidation_state = row['OS']
        
        metal_element = element(metal)
        
        if metal == 'Ru':
            sublimation = 595 / 96.485  # kJ/mol -> eV
        else:
            sublimation = metal_element.evaporation_heat / 96.485  # kJ/mol -> eV
        
        ionization = sum(metal_element.ionenergies[i] for i in range(1, oxidation_state + 1))
        
        lattice = (row['form'] - sublimation - ionization - oxygen_atomization + oxygen_electron_affinity)
        
        df.at[idx, 'sublimation_energy'] = sublimation
        df.at[idx, 'ionization_energy'] = ionization
        df.at[idx, 'lattice_energy'] = lattice
    
    return df

def plot_by_metal_row(df, save_path):
    for row in ['3d', '4d', '5d']:
        for col in ['form', 'sublimation_energy', 'ionization_energy', 'lattice_energy']:
            plt.figure(figsize=(6, 4))
            for coord in coords.index:
                zorder = coords.loc[coord, 'zorder']
                marker = coords.loc[coord, 'marker']
                color = coords.loc[coord, 'color']
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, 
                         linestyle='-', label=coord, zorder=zorder)
            
            plt.xticks(np.arange(len(metals[row])), metals[row])
            plt.xlabel("Metal Index")
            plt.ylabel(f"{col.replace('_', ' ').title()} (eV)")
            plt.legend()
            plt.tight_layout()
            png_name = f"total_bulk_{row}_{col}.png"
            plt.savefig(f"{save_path}/{png_name}")
            plt.close()
            print(f"Figure saved as {png_name}")
            
def plot_by_coordination(df, save_path):        
    for coord in coords.index:            
        for col in ['form', 'sublimation_energy', 'ionization_energy', 'lattice_energy']:
            marker = coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
            colors = cmap(np.linspace(0.0, 0.6, 3))
            plt.figure(figsize=(6, 4))
            for r, row in enumerate(['3d', '4d', '5d']):
                color = colors[r]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]           
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, 
                         linestyle='-', label=row)
            
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(f"{col.replace('_', ' ').title()} (eV)")
            plt.legend()
            plt.tight_layout()
            png_name = f"total_bulk_{coord}_{col}.png"
            plt.savefig(f"{save_path}/{png_name}")
            plt.close()
            print(f"Figure saved as {png_name}")

if __name__ == "__main__":
    # 에너지 계산
    df = calculate_energies(df)
    
    # 플롯 생성
    plot_by_metal_row(df, save_path)
    plot_by_coordination(df, save_path) 