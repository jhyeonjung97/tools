import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

print(f"\033[92m{os.getcwd()}\033[0m")

def plot_patterns_from_multiple_tsv(filenames, output, xlabel, ylabel, labels, a, b, row, fontsize):
        
    metal_rows = {
        '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
        '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
        '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
        }
    
    if row:
        indice = metal_rows[row]
        l = 9
        markers = ['>', '<', 'o', 's', 'p', 'd', 'h', '^', 'v', 'X']  # Replaced the second 'd' with 'X'
        colors = ['#d62728', '#ff7f0e', '#ffd70e', '#2ca02c', '#17becf', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
    else:
        indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metal_rows['3d'], metal_rows['4d'], metal_rows['5d'])]
        l = len(filenames)
        if '1_Tetrahedral_WZ' in os.getcwd():
            coordination = 'WZ'
            markers = ['>'] * l
            colors = plt.cm.Reds(np.linspace(0.4, 0.9, l))
        elif '2_Tetrahedral_ZB' in os.getcwd():
            coordination = 'ZB'
            markers = ['<'] * l
            colors = plt.cm.Oranges(np.linspace(0.4, 0.9, l))
        elif '3_SquarePlanar_TN' in os.getcwd():
            coordination = 'TN'
            markers = ['o'] * l
            colors = plt.cm.YlOrBr(np.linspace(0.4, 0.9, l))  # Yellow to brown gradient
        elif '4_SquarePlanar_PD' in os.getcwd():
            coordination = 'PD'
            markers = ['s'] * l
            colors = plt.cm.Greens(np.linspace(0.4, 0.9, l))
        elif '5_SquarePlanar_NB' in os.getcwd():
            coordination = 'NB'
            markers = ['p'] * l
            colors = plt.cm.Blues(np.linspace(0.4, 0.9, l))
        elif '6_Octahedral_RS' in os.getcwd():
            coordination = 'RS'
            markers = ['d'] * l
            colors = plt.cm.Purples(np.linspace(0.4, 0.9, l))
        elif '7_Pyramidal_LT' in os.getcwd():
            coordination = 'LT'
            markers = ['h'] * l  # Hexagon marker for variation
            colors = plt.cm.cool(np.linspace(0.4, 0.9, l))  # Cool colormap
        elif '8_Tetrahedral_AQ' in os.getcwd():
            coordination = 'AQ'
            markers = ['^'] * l  # Upward-pointing triangle for variation
            colors = plt.cm.summer(np.linspace(0.4, 0.9, l))  # Summer colormap
        elif '9_SquarePlanar_AU' in os.getcwd():
            coordination = 'AU'
            markers = ['v'] * l  # Downward-pointing triangle for variation
            colors = plt.cm.winter(np.linspace(0.4, 0.9, l))  # Winter colormap
        else:
            coordination = 'Unknown'
            markers = ['d'] * l  # Default to diamonds
            colors = plt.cm.Greys(np.linspace(0.4, 0.9, l))  # Default to grey
        
    merged_df = None    
    plt.figure(figsize=(a, b), dpi=300)
    
    png_filename = f"merged_{output}.png"   
    tsv_filename = f"merged_{output}.tsv"

    if l > len(labels):
        print(f"Warning: More filenames ({l}) than labels ({len(labels)}). Excess filenames will be ignored.")
        filenames = filenames[:len(labels)]

    for j, file in enumerate(filenames):
        df = pd.read_csv(file, delimiter='\t').iloc[:, 1:]
        df.columns = labels[j] if isinstance(labels[j], list) else [labels[j]]
        merged_df = pd.concat([merged_df, df], axis=1)

    # for j, column in enumerate(merged_df.columns):
    #     filtered_x = []
    #     filtered_values = []
    #     x = merged_df.index
    #     values = merged_df[column]
    #     for i, v in enumerate(values):
    #         if not np.isnan(v):
    #             filtered_x.append(i)
    #             filtered_values.append(v)
    #     if not filtered_values:
    #         print(f"No values found for pattern: {column}")
    #         continue
    #     plt.plot(filtered_x, filtered_values, marker=markers[j], color=colors[j], label=column)
        
    for j, column in enumerate(merged_df.columns):
        x = merged_df.index
        y = merged_df[column]
        plt.plot(x, y, marker=markers[j], color=colors[j], label=column)
        
    if 'hexa_ratio' in df.columns:
        plt.plot(x, [1.633]*len(x), linestyle=':', label='hexa_ratio0', color='black')
        
    if 'norm_formation' in output:
        exp_path = '/pscratch/sd/j/jiuy97/3_V_bulk/oxide/monoxides.tsv'
        if not os.path.exists(exp_path):
            exp_path = '/Users/jiuy97/Desktop/3_V_bulk/oxide/monoxides.tsv'
        if not os.path.exists(exp_path):
            exp_path = '/Users/hailey/Desktop/3_V_bulk/oxide/monoxides.tsv'
        exp_df = pd.read_csv(exp_path, delimiter='\t')
        exp_df['dH_form'] = exp_df['dH_form'] / 23.06 # 96.48
        exp_markers = {'WZ': '>', 'ZB': '<', 'TN': 'o', 
                       'PD': 's', 'NB': 'p', 'RS': 'd', 
                       'LT': 'h', 'AQ': '^', 'AU': 'v'}
        exp_colors = {'WZ': '#d62728', 'ZB': '#ff7f0e', 'TN': '#ffd70e', 
                       'PD': '#2ca02c', 'NB': '#17becf', 'RS': '#9467bd', 
                       'LT': '#8c564b', 'AQ': '#e377c2', 'AU': '#7f7f7f'}
        if row:
            for i in exp_df.index:
                if exp_df['row'][i] == row:
                    exp_marker = exp_markers.get(exp_df['Coordination'][i], '*')
                    exp_color = exp_colors.get(exp_df['Coordination'][i], '#8a8a8a')
                    plt.scatter(exp_df['numb'][i], exp_df['dH_form'][i], 
                                marker=exp_marker, color=exp_color, edgecolors=exp_color, facecolors='white', zorder=20)
        else:
            for i in exp_df.index:
                if exp_df['Coordination'][i] == coordination:
                    if exp_df['row'][i] == '3d':
                        exp_color = colors[-3]; exp_marker = markers[-3]
                    elif exp_df['row'][i] == '4d':
                        exp_color = colors[-2]; exp_marker = markers[-2]
                    elif exp_df['row'][i] == '5d':
                        exp_color = colors[-1]; exp_marker = markers[-1]
                    plt.scatter(exp_df['numb'][i], exp_df['dH_form'][i],
                                marker=exp_marker, color=exp_color, edgecolors=exp_color, facecolors='white', zorder=20)
    
    merged_df.to_csv(tsv_filename, sep='\t', float_format='%.4f')
    print(f"Merged data saved to {tsv_filename}")

    plt.xticks(np.arange(len(indice)), indice, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(prop={'size': fontsize}, ncol=1)
    plt.tight_layout()
    plt.savefig(png_filename, bbox_inches="tight")
    print(f"Figure saved as {png_filename}")
    plt.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot TSV data.')
    parser.add_argument('files', nargs='+', help='The TSV files to plot.')
    parser.add_argument('-o', '--output', type=str, default='', help="The filename for the output PNG file.")
    parser.add_argument('-x', '--xlabel', type=str, default='Element or Lattice parameter (Å)', help="xlabel")
    parser.add_argument('-y', '--ylabel', type=str, default='Energy (eV) or Charge (e)', help="ylabel")
    parser.add_argument('-l', '--labels', nargs='+', default=['Tetrahedral_WZ', 'Tetrahedral_ZB', 'SquarePlanar_TN',
                                                              'SquarePlanar_PD', 'SquarePlanar_NB', 'Octahedral_RS',
                                                              'Pyramidal_LT', 'Tetrahedral_AQ', 'SquarePlanar_AU'])
    parser.add_argument('-r', '--row', type=str, default=None)
    parser.add_argument('-a', type=float, default=8)
    parser.add_argument('-b', type=float, default=6)
    parser.add_argument('--font', type=float, default=10)
    
    args = parser.parse_args()        
    plot_patterns_from_multiple_tsv(args.files, args.output, args.xlabel, args.ylabel, args.labels, args.a, args.b, args.row, fontsize=args.font)

