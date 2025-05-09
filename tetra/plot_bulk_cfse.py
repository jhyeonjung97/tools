#!/bin/env python

import os
import socket
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

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
metal_df = pd.read_csv(f'{root}/formation_energy.csv', index_col=0)
metal_df['+3'] = metal_df['+3'] / 2
metal_df['+5'] = metal_df['+5'] / 2

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

# 컬럼 정보 설정
columns_data = [
    {'column': 'coord',    'png_name': 'coordination',        'ylabel': 'Coordination'},
    {'column': 'row',      'png_name': 'row',                 'ylabel': 'Row'},
    {'column': 'numb',     'png_name': 'number',              'ylabel': 'Number'},
    {'column': 'metal',    'png_name': 'metal',               'ylabel': 'Metal'},
    {'column': 'group',    'png_name': 'group',               'ylabel': 'Group Number'},
    {'column': 'outer_e',  'png_name': 'outer_e',             'ylabel': 'Number of Outer Electrons'},
    {'column': 'CN',       'png_name': 'coordination_number', 'ylabel': 'Coordination Number'},
    {'column': 'OS',       'png_name': 'oxidation_number',    'ylabel': 'Oxidation Number'},
    {'column': 'form',     'png_name': 'formation_energy',    'ylabel': 'Formation Energy (eV)'},
    {'column': 'volume',   'png_name': 'volume',              'ylabel': 'Volume (Å³)'},
    {'column': 'chg',      'png_name': 'bader_charge',        'ylabel': 'Bader Charge (e⁻)'},
    {'column': 'chgn',     'png_name': 'bader_charge_normalized', 'ylabel': 'normalized Bader Charge (e⁻)'},
    {'column': 'mag',      'png_name': 'magnetic_moments',    'ylabel': 'Magnetic Moments (μB)'},
    {'column': 'l_bond',   'png_name': 'bond_length',         'ylabel': 'Bond Length (Å)'},
    {'column': '-ICOHP',   'png_name': 'icohp_per_metal',     'ylabel': '-ICOHP per Metal (eV)'},
    {'column': '-ICOHPc',  'png_name': 'icohp_per_metal_c',   'ylabel': 'CN normalized -ICOHP per Metal (eV)'},
    {'column': '-ICOHPo',  'png_name': 'icohp_per_metal_o',   'ylabel': 'OS normalized -ICOHP per Metal (eV)'},
    {'column': '-ICOHPn',  'png_name': 'icohp_per_metal_n',   'ylabel': 'CN, OS normalized -ICOHP per Metal (eV)'},
    {'column': 'ICOBI',    'png_name': 'icobi_per_metal',     'ylabel': 'ICOBI per Metal'},
    {'column': 'ICOBIc',   'png_name': 'icobi_per_metal_c',   'ylabel': 'CN normalized ICOBI per Metal'},
    {'column': 'ICOBIo',   'png_name': 'icobi_per_metal_o',   'ylabel': 'OS normalized ICOBI per Metal'},
    {'column': 'ICOBIn',   'png_name': 'icobi_per_metal_n',   'ylabel': 'CN, OS normalized ICOBI per Metal'},
    {'column': '-ICOOP',   'png_name': 'icoop_per_metal',     'ylabel': '-ICOOP per Metal (eV)'},
    {'column': '-ICOOPc',  'png_name': 'icoop_per_metal_c',   'ylabel': 'CN normalized -ICOOP per Metal (eV)'},
    {'column': '-ICOOPo',  'png_name': 'icoop_per_metal_o',   'ylabel': 'OS normalized -ICOOP per Metal (eV)'},
    {'column': '-ICOOPn',  'png_name': 'icoop_per_metal_n',   'ylabel': 'CN, OS, OS normalized -ICOOP per Metal (eV)'},
    {'column': 'madelung', 'png_name': 'madelung',            'ylabel': 'Madelung Energy (Loewdin, eV)'},
    {'column': 'ion-1',    'png_name': 'ion-1',               'ylabel': 'Ionization Energy (eV)'},
    {'column': 'ion',      'png_name': 'ion',                 'ylabel': 'Ionization Energy (eV)'},
    {'column': 'ion+1',    'png_name': 'ion+1',               'ylabel': 'Ionization Energy (eV)'},
    {'column': 'ion-1c',   'png_name': 'ion_c-1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ionc',     'png_name': 'ion_c',               'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ion+1c',   'png_name': 'ion_c+1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ion-1o',   'png_name': 'ion_o-1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'iono',     'png_name': 'ion_o',               'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ion+1o',   'png_name': 'ion_o+1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ion-1n',   'png_name': 'ion_n-1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ionn',     'png_name': 'ion_n',               'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ion+1n',   'png_name': 'ion_n+1',             'ylabel': 'normalized Ionization Energy (eV)'},
    {'column': 'ionN-1',   'png_name': 'ion_accumulated-1',   'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN',     'png_name': 'ion_accumulated',     'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN+1',   'png_name': 'ion_accumulated+1',   'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN-1c',  'png_name': 'ion_accumulated_c-1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionNc',    'png_name': 'ion_accumulated_c',   'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN+1c',  'png_name': 'ion_accumulated_c+1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN-1o',  'png_name': 'ion_accumulated_o-1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionNo',    'png_name': 'ion_accumulated_o',   'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN+1o',  'png_name': 'ion_accumulated_o+1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN-1n',  'png_name': 'ion_accumulated_n-1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionNn',    'png_name': 'ion_accumulated_n',   'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'ionN+1n',  'png_name': 'ion_accumulated_n+1', 'ylabel': 'accumulated Ionization Energy (eV)'},
    {'column': 'pauling',  'png_name': 'pauling',             'ylabel': 'Pauling Electronegativity'},
    {'column': 'Natom',    'png_name': 'atomic_number',       'ylabel': 'Atomic Number'},
    {'column': 'mass',     'png_name': 'mass',                'ylabel': 'Atomic Mass (u)'},
    {'column': 'density',  'png_name': 'density',             'ylabel': 'Density (g/cm³)'},
    {'column': 'Vatom',    'png_name': 'atomic_volume',       'ylabel': 'Atomic Volume (Å³)'},
    {'column': 'dipole',   'png_name': 'dipole',              'ylabel': 'Dipole Polarizability (Å³)'},
    {'column': 'Rcoval',   'png_name': 'covalent_radius',     'ylabel': 'Covalent Radius (Å)'},
    {'column': 'Rmetal',   'png_name': 'metallic_radius',     'ylabel': 'Metallic Radius (Å)'},
    {'column': 'Rvdw',     'png_name': 'vdw_radius',          'ylabel': 'van der Waals Radius (Å)'},
    {'column': 'Tboil',    'png_name': 'boiling_point',       'ylabel': 'Boiling Point (K)'},
    {'column': 'Tmelt',    'png_name': 'melting_point',       'ylabel': 'Melting Point (K)'},
    {'column': 'Hevap',    'png_name': 'evaporation_heat',    'ylabel': 'Heat of Evaporation (kJ/mol)'},
    {'column': 'Hfus',     'png_name': 'fusion_heat',         'ylabel': 'Heat of Fusion (kJ/mol)'},
    {'column': 'Hform',    'png_name': 'formation_heat',      'ylabel': 'Heat of Formation (kJ/mol)'},
    {'column': 'cfse',     'png_name': 'cfse',                'ylabel': 'CFSE (kJ/mol)'},
    {'column': 'base_cfse', 'png_name': 'base_cfse',          'ylabel': 'Base CFSE (kJ/mol)'},
    {'column': 'jt_effect', 'png_name': 'jt_effect',          'ylabel': 'Jahn-Teller Effect (kJ/mol)'},
    {'column': 'ee_repulsion', 'png_name': 'ee_repulsion',    'ylabel': 'Electron-Electron Repulsion (kJ/mol)'},
    {'column': 'field_strength', 'png_name': 'field_strength', 'ylabel': 'Field Strength (kJ/mol)'},
    {'column': 'exchange_stabilization', 'png_name': 'exchange_stabilization', 'ylabel': 'Exchange Stabilization (kJ/mol)'},
]
columns = pd.DataFrame(columns_data).set_index('column')
columns.index.name = None

# 데이터 타입 설정
bool_cols = []
int_cols = ['CN', 'OS', 'Natom', 'n_electrons', 'd_electrons', 'outer_e', 'group']
str_cols = ['coord', 'row', 'numb', 'metal']
float_cols = [
    'form', 'volume', 'chg', 'mag', 'l_bond', 'madelung', 
    '-ICOHP', 'ICOBI', '-ICOOP', '-ICOHPo', 'ICOBIo', '-ICOOPo', '-ICOHPn', 'ICOBIc', '-ICOOPc', 
    'ion-1', 'ion', 'ion+1', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1', 'pauling',
    'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
    'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform', 
    'base_cfse', 'ee_repulsion', 'jt_effect', 'field_strength', 'cfse', 'exchange_stabilization',
    ]
    		
def plot_by_metal_row(df, save_path):
    for row in ['3d', '4d', '5d']:
        for col in columns.index:
            if col in str_cols or col in bool_cols:
                continue
            plt.figure(figsize=(8, 6))
            for coord in coords.index:
                zorder = coords.loc[coord, 'zorder']
                marker = coords.loc[coord, 'marker']
                color = coords.loc[coord, 'color']
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, 
                         linestyle='-', label=coord, zorder=zorder)
                if col == 'form':
                    for m, metal in enumerate(metals[row]):
                        plt.scatter(m, metal_df.loc[metal, coord]/23.06, marker=marker, edgecolors=color, facecolors='white', zorder=zorder)
                    
                    # d0-d5, d5-d10을 잇는 선 추가
                    if coord in ['+3', '+4', '+5', '+6']:
                        # +3 ~ +6 coord는 기존 방식대로
                        d0_mask = (subset['n_electrons'] == 0)
                        d5_mask = (subset['n_electrons'] == 5)
                        d10_mask = (subset['n_electrons'] == 10)
                        if d0_mask.any() and d5_mask.any():
                            d0_idx = subset[d0_mask].index[0]
                            d5_idx = subset[d5_mask].index[0]
                            plt.plot([subset.loc[d0_idx, 'numb'], subset.loc[d5_idx, 'numb']], 
                                    [subset.loc[d0_idx, 'form'], subset.loc[d5_idx, 'form']],
                                    '--', color=color, zorder=zorder-1, linewidth=0.5, dashes=(12, 5))
                        if d5_mask.any() and d10_mask.any():
                            d5_idx = subset[d5_mask].index[0]
                            d10_idx = subset[d10_mask].index[0]
                            plt.plot([subset.loc[d5_idx, 'numb'], subset.loc[d10_idx, 'numb']], 
                                    [subset.loc[d5_idx, 'form'], subset.loc[d10_idx, 'form']],
                                    '--', color=color, zorder=zorder-1, linewidth=0.5, dashes=(12, 5))
                    else:
                        # WZ ~ LT coord는 최소값으로 연결
                        tetra_coords = ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']
                        tetra_subset = df[(df['coord'].isin(tetra_coords)) & (df['row'] == row)]
                        
                        # d0, d5, d10에 해당하는 최소 formation energy 찾기
                        d0_mask = tetra_subset['n_electrons'] == 0
                        d5_mask = tetra_subset['n_electrons'] == 5
                        d10_mask = tetra_subset['n_electrons'] == 10
                        
                        if d0_mask.any() and d5_mask.any():
                            d0_min = tetra_subset[d0_mask]['form'].min()
                            d5_min = tetra_subset[d5_mask]['form'].min()
                            
                            # d0, d5에 해당하는 metal 찾기
                            d0_metal = tetra_subset[(tetra_subset['form'] == d0_min) & d0_mask]['metal'].iloc[0]
                            d5_metal = tetra_subset[(tetra_subset['form'] == d5_min) & d5_mask]['metal'].iloc[0]
                            
                            # 선 그리기
                            d0_idx = metals[row].index(d0_metal)
                            d5_idx = metals[row].index(d5_metal)
                            
                            plt.plot([d0_idx, d5_idx], [d0_min, d5_min],
                                    '--', color='red', zorder=zorder-1, linewidth=0.5, dashes=(12, 5))
                        
                        if d5_mask.any() and d10_mask.any():
                            d5_min = tetra_subset[d5_mask]['form'].min()
                            d10_min = tetra_subset[d10_mask]['form'].min()
                            
                            # d5, d10에 해당하는 metal 찾기
                            d5_metal = tetra_subset[(tetra_subset['form'] == d5_min) & d5_mask]['metal'].iloc[0]
                            d10_metal = tetra_subset[(tetra_subset['form'] == d10_min) & d10_mask]['metal'].iloc[0]
                            
                            # 선 그리기
                            d5_idx = metals[row].index(d5_metal)
                            d10_idx = metals[row].index(d10_metal)
                            
                            plt.plot([d5_idx, d10_idx], [d5_min, d10_min],
                                    '--', color='red', zorder=zorder-1, linewidth=0.5, dashes=(12, 5))
                            
            plt.xticks(np.arange(len(metals[row])), metals[row])
            plt.xlabel("Metal Index")
            plt.ylabel(columns.loc[col, 'ylabel'])
            plt.xlim(-0.5, 12.5)  # x축 범위를 0부터 13까지로 설정
            plt.legend()
            plt.tight_layout()
            png_name = f"total_bulk_{row}_{columns.loc[col, 'png_name']}.png"
            plt.savefig(f"{save_path}/{png_name}", dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")
            
def plot_by_coordination(df, save_path):        
    for coord in coords.index:            
        for col in columns.index:
            if col in str_cols or col in bool_cols:
                continue
            marker = coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
            colors = cmap(np.linspace(0.0, 0.6, 3))
            plt.figure(figsize=(8, 6))
            for r, row in enumerate(['3d', '4d', '5d']):
                color = colors[r]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]           
                plt.plot(subset['numb'], subset[col], marker=marker, color=color, 
                         linestyle='-', label=row)
                if col == 'form':
                    for m, metal in enumerate(metals[row]):
                        plt.scatter(m, metal_df.loc[metal, coord]/23.06, marker=marker, edgecolors=color, facecolors='white')
                    
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(columns.loc[col, 'ylabel'])
            plt.xlim(-0.5, 12.5)  # x축 범위를 0부터 13까지로 설정
            plt.legend()
            plt.tight_layout()
            png_name = f"total_bulk_{coord}_{columns.loc[col, 'png_name']}.png"
            plt.savefig(f"{save_path}/{png_name}", dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")

if __name__ == "__main__":
    plot_by_metal_row(df, save_path)
    plot_by_coordination(df, save_path) 

cfse_ylabels = {
    'd_electrons': 'Number of d-electrons',
    'outer_e': 'Number of Outer Electrons',
    'group': 'Group Number',
    'base_cfse': 'Base CFSE (eV)',
    'ee_repulsion': 'Electron-Electron Repulsion (eV)',
    'jt_effect': 'Jahn-Teller Effect (eV)',
    'exchange_stabilization': 'Exchange Stabilization (eV)',
    'field_strength': 'Crystal Field Strength',
    'cfse': 'Total CFSE (eV)'
} 