import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

# Define the rows and spins
rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
spins = {'LS': '#ff7f0e', 'IS': '#279ff2', 'HS': '#9467bd'}
ms_spins ={'MS(LS)': '#ff7f0e', 'MS(IS)': '#279ff2', 'MS(HS)': '#9467bd'}
dzs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

E_H2O = -14.23919983
E_H2 = -6.77409008
E_OH = E_H2O - E_H2/2
E_O = E_H2O - E_H2

# Energies and corrections
nitrogen_E = -16.64503942 # eV, DFT
nitrogen_TS = 0.635139 # eV, at 298.15 K, 1 atm
nitrogen_ZPE = 0.096279 # eV, at 298.15 K, 1 atm
nitrogen = (nitrogen_E - nitrogen_TS + nitrogen_ZPE) / 2

carbon = -9.357363435 # eV, DFT

metal_path = '/pscratch/sd/j/jiuy97/6_MNC/gas/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)

def main():
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals):
            df = pd.DataFrame()
            # df_rel = pd.DataFrame()
            df_mag = pd.DataFrame()
            df_relaxed = pd.DataFrame()
            # df_relaxed_rel = pd.DataFrame()
            df_relaxed_mag = pd.DataFrame()

            df_O = pd.DataFrame()
            # df_O_rel = pd.DataFrame()
            df_O_mag = pd.DataFrame()
            df_O_relaxed = pd.DataFrame()
            # df_O_relaxed_rel = pd.DataFrame()
            df_O_relaxed_mag = pd.DataFrame()

            df_OH = pd.DataFrame()
            # df_OH_rel = pd.DataFrame()
            df_OH_mag = pd.DataFrame()
            df_OH_relaxed = pd.DataFrame()
            # df_OH_relaxed_rel = pd.DataFrame()
            df_OH_relaxed_mag = pd.DataFrame()

            tsv_filename = f'{row_key}_{m+2}{metal}_Ef.tsv'
            png_filename = f'{row_key}_{m+2}{metal}_Ef.png'
            # tsv_rel_filename = f'{row_key}_{m+2}{metal}_rel.tsv'
            # png_rel_filename = f'{row_key}_{m+2}{metal}_rel.png'
            tsv_mag_filename = f'{row_key}_{m+2}{metal}_mag.tsv'
            png_mag_filename = f'{row_key}_{m+2}{metal}_mag.png'

            tsv_O_filename = f'{row_key}_{m+2}{metal}_ads_O.tsv'
            png_O_filename = f'{row_key}_{m+2}{metal}_ads_O.png'
            # tsv_O_rel_filename = f'{row_key}_{m+2}{metal}_rel_O.tsv'
            # png_O_rel_filename = f'{row_key}_{m+2}{metal}_rel_O.png'
            tsv_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.tsv'
            png_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.png'

            tsv_OH_filename = f'{row_key}_{m+2}{metal}_OHads.tsv'
            png_OH_filename = f'{row_key}_{m+2}{metal}_OHads.png'
            # tsv_OH_rel_filename = f'{row_key}_{m+2}{metal}_rel_OH.tsv'
            # png_OH_rel_filename = f'{row_key}_{m+2}{metal}_rel_OH.png'
            tsv_OH_mag_filename = f'{row_key}_{m+2}{metal}_mag_OH.tsv'
            png_OH_mag_filename = f'{row_key}_{m+2}{metal}_mag_OH.png'

            for spin in spins.keys():
                path_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/*_{metal}/*_{spin}'
                matching_paths = glob.glob(path_pattern)
                path_O_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/1_O/*_{metal}/*_{spin}'
                matching_O_paths = glob.glob(path_O_pattern)
                path_OH_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/1_O/*_{metal}/*_{spin}'
                matching_OH_paths = glob.glob(path_OH_pattern)
                
                path = None
                path_O = None
                path_OH = None
                
                for i, dz in enumerate(dzs):
                    for path in matching_paths:
                        atoms_path = os.path.join(path, f'{i}_', 'final_with_calculator.json')
                        if os.path.exists(atoms_path):
                            atoms = read(atoms_path)
                            energy = atoms.get_total_energy()
                            formation_energy = energy - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
                            df.at[dz, spin] = formation_energy
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_mag.at[dz, spin] = magmoms[atom.index]
                            except:
                                df_mag.at[dz, spin] = 0
                        else:
                            df.at[dz, spin] = np.nan
                            df_mag.at[dz, spin] = np.nan

                    for path_O in matching_O_paths:
                        atoms_path = os.path.join(path_O, f'{i}_', 'final_with_calculator.json')
                        if os.path.exists(atoms_path): # and energy:
                            atoms = read(atoms_path)
                            energy_O = atoms.get_total_energy()
                            adsorption_energy = energy_O - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen - E_O
                            # adsorption_energy = energy_O - energy - E_O
                            df_O.at[dz, spin] = adsorption_energy
                            energy = None
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_O_mag.at[dz, spin] = magmoms[atom.index]
                            except:
                                df_O_mag.at[dz, spin] = 0
                        else:
                            df_O.at[dz, spin] = np.nan
                            df_O_mag.at[dz, spin] = np.nan
                            
                    for path_OH in matching_OH_paths:
                        atoms_path = os.path.join(path_OH, f'{i}_', 'final_with_calculator.json')
                        if os.path.exists(atoms_path): # and energy:
                            atoms = read(atoms_path)
                            energy_OH = atoms.get_total_energy()
                            adsorption_energy = energy_OH - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen - E_OH
                            # adsorption_energy = energy_O - energy - E_O
                            df_OH.at[dz, spin] = adsorption_energy
                            energy = None
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_OH_mag.at[dz, spin] = magmoms[atom.index]
                            except:
                                df_OH_mag.at[dz, spin] = 0
                        else:
                            df_OH.at[dz, spin] = np.nan
                            df_OH_mag.at[dz, spin] = np.nan
                            
                if path:
                    relaxed_path = os.path.join(path, 'relaxed', 'final_with_calculator.json')
                    if os.path.exists(relaxed_path):
                        atoms = read(relaxed_path)
                        zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                        zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                        dz_relaxed = abs(zN - zM)
                        energy = atoms.get_total_energy()
                        formation_energy = energy - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
                        df_relaxed.at[dz_relaxed, spin] = formation_energy
                        try:
                            magmoms = atoms.get_magnetic_moments()
                            for atom in atoms:
                                if atom.symbol not in ['N', 'C', 'O', 'H']:
                                    df_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                        except:
                            df_relaxed_mag.at[dz_relaxed, spin] = 0
                        
                if path_O:
                    relaxed_O_path = os.path.join(path_O, 'relaxed', 'final_with_calculator.json')
                    if os.path.exists(relaxed_O_path): # and energy:
                        atoms = read(relaxed_O_path)
                        zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                        zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                        dz_relaxed = abs(zN - zM)
                        energy_O = atoms.get_total_energy()
                        adsorption_energy = energy_O - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen - E_O
                        # adsorption_energy = energy_O - energy - E_O
                        df_O_relaxed.at[dz_relaxed, spin] = adsorption_energy
                        energy = None
                        try:
                            magmoms = atoms.get_magnetic_moments()
                            for atom in atoms:
                                if atom.symbol not in ['N', 'C', 'O', 'H']:
                                    df_O_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                        except:
                            df_O_relaxed_mag.at[dz_relaxed, spin] = 0
                            
                if path_OH:
                    relaxed_OH_path = os.path.join(path_OH, 'relaxed', 'final_with_calculator.json')
                    if os.path.exists(relaxed_OH_path): # and energy:
                        atoms = read(relaxed_OH_path)
                        zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                        zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                        dz_relaxed = abs(zN - zM)
                        energy_OH = atoms.get_total_energy()
                        adsorption_energy = energy_OH - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen - E_OH
                        # adsorption_energy = energy_O - energy - E_O
                        df_OH_relaxed.at[dz_relaxed, spin] = adsorption_energy
                        energy = None
                        try:
                            magmoms = atoms.get_magnetic_moments()
                            for atom in atoms:
                                if atom.symbol not in ['N', 'C', 'O', 'H']:
                                    df_OH_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                        except:
                            df_OH_relaxed_mag.at[dz_relaxed, spin] = 0
            
            path_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/*_{metal}'
            matching_paths = glob.glob(path_pattern)

            for path in matching_paths:
                spin_tsv = os.path.join(path, 'lowest.tsv')
                if os.path.exists(spin_tsv):
                    spin_df = pd.read_csv(spin_tsv, sep='\t')
                    spin_df.set_index('dz', inplace=True)
                    for i, dz in enumerate(dzs):
                        if len(spin_df) > 0:
                            ms = spin_df.loc[dz, 'spin_state']
                        else:
                            continue

                        atoms_path = os.path.join(path, 'most_stable', f'{i}_', 'final_with_calculator.json')
                        if os.path.exists(atoms_path):
                            atoms = read(atoms_path)
                            energy = atoms.get_total_energy()
                            formation_energy = energy - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
                            df.at[dz, f'MS({ms})'] = formation_energy
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_mag.at[dz, f'MS({ms})'] = magmoms[atom.index]
                            except:
                                df_mag.at[dz, f'MS({ms})'] = 0
                        else:
                            df.at[dz, f'MS({ms})'] = np.nan
                            df_mag.at[dz, f'MS({ms})'] = np.nan
                            
            # relative(df, df_rel)
            # relative(df_O, df_O_rel)
            # relative(df_OH, df_OH_rel)
            # relative(df_relaxed, df_relaxed_rel)
            # relative(df_O_relaxed, df_O_relaxed_rel)
            # relative(df_OH_relaxed, df_OH_relaxed_rel)

            combining(df=df, df_relaxed=df_relaxed, tsv_filename=tsv_filename)
            # combining(df=df_rel, df_relaxed=df_relaxed_rel, tsv_filename=tsv_rel_filename)
            combining(df=df_mag, df_relaxed=df_relaxed_mag, tsv_filename=tsv_mag_filename)
            
            plotting(df=df, df_relaxed=df_relaxed, dzs=dzs, spins=spins, 
                     ylabel='Formation energy (eV)', png_filename=png_filename)
            # plotting(df=df_rel, df_relaxed=df_relaxed_rel, dzs=dzs, spins=spins, color='black', 
            #          ylabel='Spin crossover energy (eV)', png_filename=png_rel_filename)
            plotting(df=df_mag, df_relaxed=df_relaxed_mag, dzs=dzs, spins=spins, 
                     ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                     ylabel='Magnetic Moments (uB)', png_filename=png_mag_filename)
            
            if path_O: 
                combining(df=df_O, df_relaxed=df_O_relaxed, tsv_filename=tsv_O_filename)
                # combining(df=df_O_rel, df_relaxed=df_O_relaxed_rel, tsv_filename=tsv_O_rel_filename)
                combining(df=df_O_mag, df_relaxed=df_O_relaxed_mag, tsv_filename=tsv_O_mag_filename)
                
                plotting(df=df_O, df_relaxed=df_O_relaxed, dzs=dzs, spins=spins, 
                         ylabel='Formation energy (eV)', png_filename=png_O_filename)
                # plotting(df=df_O_rel, df_relaxed=df_O_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                #          ylabel='Spin crossover energy (eV)', png_filename=png_O_rel_filename)
                plotting(df=df_O_mag, df_relaxed=df_O_relaxed_mag, dzs=dzs, spins=spins, 
                         ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                         ylabel='Magnetic Moments (uB)', png_filename=png_O_mag_filename)
                
            if path_OH: 
                combining(df=df_OH, df_relaxed=df_OH_relaxed, tsv_filename=tsv_OH_filename)
                # combining(df=df_OH_rel, df_relaxed=df_OH_relaxed_rel, tsv_filename=tsv_OH_rel_filename)
                combining(df=df_OH_mag, df_relaxed=df_OH_relaxed_mag, tsv_filename=tsv_OH_mag_filename)
                
                plotting(df=df_OH, df_relaxed=df_OH_relaxed, dzs=dzs, spins=spins, 
                         ylabel='Formation energy (eV)', png_filename=png_OH_filename)
                # plotting(df=df_OH_rel, df_relaxed=df_OH_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                #          ylabel='Spin crossover energy (eV)', png_filename=png_OH_rel_filename)
                plotting(df=df_OH_mag, df_relaxed=df_OH_relaxed_mag, dzs=dzs, spins=spins, 
                         ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                         ylabel='Magnetic Moments', png_filename=png_OH_mag_filename)
                
# def relative(df, df_rel):
#     if 'HS' in df.columns and 'LS' in df.columns:
#         df_rel['HS-LS'] = df['HS'] - df['LS']

def combining(df, df_relaxed, tsv_filename):
    combined_df = pd.concat([df, df_relaxed])
    combined_df.to_csv(tsv_filename, sep='\t', float_format='%.2f')
    print(f"Data saved to {tsv_filename}")

def plot_smooth_line(x, y, color):
    try:
        x_new = np.linspace(min(x), max(x), 300)
        if len(x) > 3:
            spl = make_interp_spline(x, y, k=3)  # Smoothing spline
        else:
            spl = make_interp_spline(x, y, k=2)
        y_smooth = spl(x_new)
        plt.plot(x_new, y_smooth, color=color, alpha=0.3, zorder=1)
        plt.scatter(x, y, marker='s', edgecolors=color, facecolors='white', alpha=0.3, zorder=2)
        return y_smooth
    except ValueError as e:
        print(f"Error while creating spline: {e}")
        return None

def plotting(df, df_relaxed, dzs, spins, ylabel, png_filename, ymin=None, ymax=None, yticks=None, color=None):
    plt.figure(figsize=(4, 3))
    df_smooth_y = pd.DataFrame()
    for column in df.columns:
        filtered_df = df[column].dropna()
        if not filtered_df.empty:
            x = filtered_df.index
            y = filtered_df.values
            if 'MS' in column:
                plt.scatter(x, y, marker='x', color=ms_spins.get(column, 'black'), zorder=3)
            else:
                df_smooth_y[column] = plot_smooth_line(x, y, color or spins.get(column, 'black'))
    
    min_values = df_smooth_y.min(axis=1).to_numpy()
    min_columns = df_smooth_y.idxmin(axis=1).to_numpy()
    x_new = np.linspace(0.0, 1.2, 300)
    
    start_idx = 0
    current_column = min_columns[0]
    for i in range(1, len(min_columns)):
        if min_columns[i] != current_column:
            plt.plot(x_new[start_idx:i], min_values[start_idx:i], color=spins.get(current_column, 'black'), zorder=5)
            start_idx = i
            current_column = min_columns[i]
    plt.plot(x_new[start_idx:], min_values[start_idx:], color=spins.get(current_column, 'black'), zorder=5)
    
    for column in df_relaxed.columns:
        filtered_df = df_relaxed[column].dropna()
        if not filtered_df.empty:
            x = filtered_df.index
            y = filtered_df.values
            plt.scatter(x, y, marker='s', color=color or spins.get(column, 'black'), alpha=0.3, zorder=4)
    if color:
        plt.axhline(y=0.0, color='blue', linestyle='--', zorder=0)
        plt.axhline(y=0.8, color='red', linestyle='--', zorder=0)
    plt.xticks(dzs)
    plt.xlabel('dz (Å)')
    plt.ylabel(ylabel)
    if ymin and ymax:
        plt.ylim(ymin, ymax)
    if yticks is not None:
        plt.yticks(yticks)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # Fix to 0.0 format
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # Fix to 0.0 format
    # plt.legend(labelspacing=0.3)
    plt.tight_layout()
    plt.savefig(png_filename, bbox_inches="tight")
    print(f"Figure saved as {png_filename}")
    plt.close()

if __name__ == '__main__':
    main()
    custom_legend = [
        Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#ff7f0e', 
               label='LS (fixed, w/ nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#279ff2', 
               label='IS (fixed, w/ nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#9467bd', 
               label='HS (fixed, w/ nupdown)', markersize=8, linestyle='None'),
    
        Line2D([0], [0], marker='s', color='#ff7f0e', 
               label='LS (relaxed, w/ nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='s', color='#279ff2', 
               label='IS (relaxed, w/ nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='s', color='#9467bd', 
               label='HS (relaxed, w/ nupdown)', markersize=8, linestyle='None'),
    
        Line2D([0], [0], marker='x', color='#ff7f0e', 
               label='LS (fixed, w/o nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='x', color='#279ff2', 
               label='IS (fixed, w/o nupdown)', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='x', color='#9467bd', 
               label='HS (fixed, w/o nupdown)', markersize=8, linestyle='None'),
    ]
    fig, ax = plt.subplots()
    ax.legend(handles=custom_legend)
    ax.axis('off')
    png_filename = "custom_legend_markers_only.png"  # Update with your file path
    plt.savefig(png_filename, bbox_inches="tight")
    print(f"Figure saved as {png_filename}")
    plt.close()