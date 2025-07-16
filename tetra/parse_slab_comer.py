#!/bin/env python

import os
import socket
import numpy as np
import pandas as pd
from ase.io import read, write
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 서버 주소 가져오기
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/8_V_slab'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/8_V_slab'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/8_V_slab'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")
save_path = os.path.join(root, 'figures')
os.makedirs(save_path, exist_ok=True)

# 참조 에너지 값들
h2o = -14.23919983
h2 = -6.77409008

zpeh2o = 0.558
zpeh2 = 0.273

cph2o = 0.10
cph2 = 0.09

gh2 = h2 + zpeh2 + cph2
gh2o = h2o + zpeh2o + cph2o
go = gh2o - gh2
goh = gh2o - gh2/2

coords_data = [
    {'coord': '+3', 'CN': 6, 'OS': 3, 'MN': 4, 'coord_dir': '1_Octahedral_+3_012', 'zorder': 1, 'marker': 'o', 'color': 'darkorange'},
    {'coord': '+4', 'CN': 6, 'OS': 4, 'MN': 2, 'coord_dir': '2_Octahedral_+4_100', 'zorder': 2, 'marker': 's', 'color': 'gold'},
    {'coord': '+4', 'CN': 6, 'OS': 4, 'MN': 2, 'coord_dir': '3_Octahedral_+4_110', 'zorder': 2, 'marker': 's', 'color': 'silver'},
    {'coord': '+5', 'CN': 6, 'OS': 5, 'MN': 4, 'coord_dir': '4_Octahedral_+5_100', 'zorder': 3, 'marker': '^', 'color': 'dodgerblue'},
    {'coord': '+6', 'CN': 6, 'OS': 6, 'MN': 1, 'coord_dir': '5_Octahedral_+6_001', 'zorder': 4, 'marker': 'v', 'color': 'deepskyblue'},
]

coords = pd.DataFrame(coords_data).set_index('coord_dir')
coords.index.name = None

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
}

def main():
    df = pd.DataFrame(columns=['coord', 'row', 'numb', 'metal', 'surfox', 'o_energy', 'oh_energy'])
    
    dir2surface = {
    '1_Octahedral_+3_012': 2.833,
    '2_Octahedral_+4_100': 4.000,
    '3_Octahedral_+4_110': 3.333,
    '4_Octahedral_+5_100': 5.000,
    '5_Octahedral_+6_001': 6.000,
    }

    for coord_dir in coords.index:
        coord = coords.loc[coord_dir, 'coord']
        surfox = dir2surface[coord_dir]
        for row in ['3d', '4d', '5d']:
            metals_list = metals[row]
            for m, metal in enumerate(metals_list):
                m_str = str(m).zfill(2)
                item = coord + row + m_str
                path = os.path.join(root, 'comer', coord_dir, row, f'{m_str}_{metal}')
                clean_path = os.path.join(path, 'empty_slab.json')
                o_path = os.path.join(path, 'O.json')
                oh_path = os.path.join(path, 'HO.json')
                clean_energy = get_energy(clean_path)
                o_energy = get_energy(o_path)
                oh_energy = get_energy(oh_path)
                df.loc[item, ['coord', 'row', 'numb', 'metal', 'surfox']] = coord, row, m_str, metal, surfox
                if clean_energy is not None and o_energy is not None:
                    o_ads_energy = o_energy - clean_energy - go + 0.28931784
                    df.loc[item, 'o_energy'] = o_ads_energy
                if clean_energy is not None and oh_energy is not None:
                    oh_ads_energy = oh_energy - clean_energy - goh + 0.46951875
                    df.loc[item, 'oh_energy'] = oh_ads_energy

    # 결과 저장
    df = df.dropna(subset=['o_energy', 'oh_energy'], how='all')
    df.to_csv(f'{save_path}/slab_data_comer.csv', sep=',')
    df[['o_energy', 'oh_energy']] = df[['o_energy', 'oh_energy']].astype(float).round(2)
    df.to_csv(f'{save_path}/slab_data_comer.tsv', sep='\t', float_format='%.2f')
    
    # coord, row, numb 순으로 정렬 (numb를 int로 변환)
    df['numb'] = df['numb'].astype(int)
    df = df.sort_values(['coord', 'row', 'numb'])
    print(df)
    
    # 그래프 그리기
    # plot_by_metal_row(df, save_path)
    # plot_by_coordination(df, save_path)

def get_energy(atoms_path):
    if os.path.exists(atoms_path):
        atoms = read(atoms_path)
        return atoms.get_total_energy()
    return None

def plot_by_metal_row(df, save_path):
    """각 금속 행(row)별로 흡착 에너지를 그래프로 그립니다."""
    for row in ['fm', '3d', '4d', '5d']:
        metals_list = metals[row]
        for ads in ['o_energy', 'oh_energy']:
            plt.figure(figsize=(4, 3))
            for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT', 'WW', 'ZZ']:
                zorder = coords.loc[coord, 'zorder']
                marker = coords.loc[coord, 'marker']
                color = coords.loc[coord, 'color']
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                # metals_list 순서에 맞게 reindex
                subset = subset.set_index('metal').reindex(metals_list).reset_index()
                plt.plot(np.arange(len(metals_list)), subset[ads], marker=marker, color=color, 
                        linestyle='-', label=coord, zorder=zorder)
            plt.xticks(np.arange(len(metals_list)), metals_list)
            plt.xlabel("Metal Index")
            plt.ylabel(f"{ads.split('_')[0].upper()} Adsorption Energy (eV)")
            plt.xlim(-0.5, 12.5)
            plt.legend()
            plt.tight_layout()
            png_name = f"slab_comer_{row}_{ads}.png"
            plt.savefig(f"{save_path}/{png_name}", transparent=True, dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")

def plot_by_coordination(df, save_path):
    """각 좌표계(coord)별로 흡착 에너지를 그래프로 그립니다."""
    for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT', 'WW', 'ZZ']:
        for ads in ['o_energy', 'oh_energy']:
            marker = coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
            colors = cmap(np.linspace(0.0, 0.6, 3))
            plt.figure(figsize=(4, 3))
            for r, row in enumerate(['3d', '4d', '5d']):
                metals_list = metals[row]
                color = 'lightgray' if row == 'fm' else colors[r]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                subset = subset.set_index('metal').reindex(metals_list).reset_index()
                plt.plot(np.arange(len(metals_list)), subset[ads], marker=marker, color=color, 
                        linestyle='-', label=row)
            plt.xticks(np.arange(len(metals['3d'])), [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])])
            plt.xlabel("Metal Index")
            plt.ylabel(f"{ads.split('_')[0].upper()} Adsorption Energy (eV)")
            plt.xlim(-0.5, 12.5)
            plt.legend()
            plt.tight_layout()
            png_name = f"slab_comer_{coord}_{ads}.png"
            plt.savefig(f"{save_path}/{png_name}", transparent=True, dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")

if __name__ == '__main__':
    main() 