import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

coords = ['WZ', 'ZB', 'LT', 'TN', 'NB', 'RS']
coord_dirs = ['1_Tetrahedral_WZ', '2_Tetrahedral_ZB', '3_Pyramidal_LT',
              '4_Square_Planar_TN', '5_Square_Planar_33', '6_Octahedral_RS']
colors = ['#d62728', '#ff7f0e', '#ffd70e', '#2ca02c', '#279ff2', '#9467bd']
color_ranges = [plt.cm.Reds(np.linspace(0.3, 0.9, 3)),
                plt.cm.Oranges(np.linspace(0.3, 0.9, 3)),
                plt.cm.Wistia(np.linspace(0.3, 0.9, 3)),
                plt.cm.Greens(np.linspace(0.3, 0.9, 3)),
                plt.cm.Blues(np.linspace(0.3, 0.9, 3)),
                plt.cm.Purples(np.linspace(0.3, 0.9, 3))]
markers = ['>', '<', 'o', 's', 'p', 'd']
stochiometries = [8, 8, 8, 8, 12, 8]

rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}
row_dirs = ['1_afm', '4d', '5d']

indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(rows['3d'], rows['4d'], rows['5d'])]

bulk_path = '/pscratch/sd/j/jiuy97/3_V_bulk'
slab_path = '/pscratch/sd/j/jiuy97/4_V_slab'

for i in range(6):
    coord = coords[i]
    coord_dir = coord_dirs[i]
    color = colors[i]
    color_range = color_ranges[i]
    marker = markers[i]
    stochiometry = stochiometries[i]
    combined_df = pd.DataFrame()
    
    for j in range(3):
        row_key = list(rows.keys())[j]
        row = rows[row_key]
        row_dir = row_dirs[j]
        
        bulk_e_path = os.path.join(bulk_path, coord_dir, row_dir, 'energy_norm_energy.tsv')
        slab_e_path = os.path.join(slab_path, coord_dir, row_key, 'energy_energy.tsv')
        area_e_path = os.path.join(slab_path, coord_dir, row_key, 'energy_area.tsv')
        
        if os.path.exists(bulk_e_path) and os.path.exists(slab_e_path) and os.path.exists(area_e_path):
            bulk_df = pd.read_csv(bulk_e_path, delimiter='\t').iloc[:, 1:]
            slab_df = pd.read_csv(slab_e_path, delimiter='\t').iloc[:, 1:]
            area_df = pd.read_csv(area_e_path, delimiter='\t').iloc[:, 1:]
            
            surface_df = pd.DataFrame(index=bulk_df.index, columns=bulk_df.columns)
            
            for k in range(len(bulk_df)):
                if not (pd.isna(slab_df.iloc[k, 0]) or pd.isna(bulk_df.iloc[k, 0]) or pd.isna(area_df.iloc[k, 0])):
                    if coord == 'NB' and row_key == '3d':
                        surface_df.iloc[k, 0] = (slab_df.iloc[k, 0] - 24 * bulk_df.iloc[k, 0]) / (2 * area_df.iloc[k, 0])
                    # elif coord == 'WZ' and row_key == '3d':                        
                    else:
                        surface_df.iloc[k, 0] = (slab_df.iloc[k, 0] - stochiometry * bulk_df.iloc[k, 0]) / (2 * area_df.iloc[k, 0])
                else:
                    surface_df.iloc[k, 0] = np.nan
                    
            surface_df['energy'] = surface_df['energy'].astype(float)
            
            max_value = surface_df.max().max()
            min_value = surface_df.min().min()
            # print(f"Max, Min for {coord}_{row_key}: {max_value}, {min_value}")
            print("Max, Min for {}_{}: {:.4f}, {:.4f}".format(coord, row_key, max_value, min_value))

            combined_df = pd.concat([combined_df, surface_df.rename(columns={'energy': f'{coord}_{row_key}'})], axis=1)
            
            png_filename = f"surface_{i}{coord}_{row_key}.png"
            tsv_filename = f"surface_{i}{coord}_{row_key}.tsv"

            plt.figure(figsize=(8, 6))
            x = range(len(surface_df['energy']))
            plt.axhline(y=0, color='gray', linestyle='--')
            plt.plot(x, surface_df['energy'], marker=marker, color=color, label=f'{coord}_{row_key}')
            surface_df.to_csv(tsv_filename, sep='\t', float_format='%.4f')
            print(f"Merged data saved to {tsv_filename}")
            plt.xticks(np.arange(len(row)), row)
            plt.xlabel('Metal (MO)')
            plt.ylabel('Surface energy (eV/A^2)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(png_filename, bbox_inches="tight")
            print(f"Figure saved as {png_filename}")
            plt.close()

    png_filename_combined = f"surface_{i}{coord}.png"
    tsv_filename_combined = f"surface_{i}{coord}.tsv"
    
    plt.figure(figsize=(8, 6))
    for m, column in enumerate(combined_df.columns):
        x = range(len(combined_df[column]))
        plt.axhline(y=0, color='gray', linestyle='--')
        plt.plot(x, combined_df[column], marker=marker, color=color_range[m], label=column)
    combined_df.to_csv(tsv_filename_combined, sep='\t', float_format='%.4f')
    print(f"Combined data saved to {tsv_filename_combined}")
    plt.xticks(np.arange(len(indice)), indice)
    plt.xlabel('Metal (MO)')
    plt.ylabel('Surface energy (eV/A^2)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(png_filename_combined, bbox_inches="tight")
    print(f"Figure saved as {png_filename_combined}")
    plt.close()