import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np

metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
spins = {'LS': '#ff7f0e', 'IS': '#279ff2', 'HS': '#9467bd'}
dzs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
path = '/pscratch/sd/j/jiuy97/6_MNC/figure'

def main():
    for m, metal in enumerate(metals):
        path_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/figure/*d_*{metal}_Ef.tsv'
        matching_O_paths = glob.glob(path_O_pattern)

        for path_O in matching_O_paths:
            print(path_O)
            for i, dz in enumerate(dzs):
                atoms_O_path = os.path.join(path, f'{i}_', 'final_with_calculator.json')
                if os.path.exists(atoms_O_path):
                        
        df = pd.DataFrame()
        df_mag = pd.DataFrame()
        df_relaxed = pd.DataFrame()
        df_relaxed_mag = pd.DataFrame()

        df_O = pd.DataFrame()
        df_O_mag = pd.DataFrame()
        df_O_relaxed = pd.DataFrame()
        df_O_relaxed_mag = pd.DataFrame()

        tsv_filename = f'{row_key}_{m+2}{metal}_Ef.tsv'
        png_filename = f'{row_key}_{m+2}{metal}_Ef.png'
        tsv_mag_filename = f'{row_key}_{m+2}{metal}_mag.tsv'
        png_mag_filename = f'{row_key}_{m+2}{metal}_mag.png'

        tsv_O_filename = f'{row_key}_{m+2}{metal}_Oads.tsv'
        png_O_filename = f'{row_key}_{m+2}{metal}_Oads.png'
        tsv_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.tsv'
        png_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.png'

        for spin in spins.keys():
            path_O_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/1_O/*_{metal}/*_{spin}'
            matching_O_paths = glob.glob(path_O_pattern)

            for path_O in matching_O_paths:
                print(path_O)
                for i, dz in enumerate(dzs):
                    atoms_O_path = os.path.join(path, f'{i}_', 'final_with_calculator.json')
                    if os.path.exists(atoms_O_path):
                        atoms_O = read(atoms_O_path)
                        energy_O = atoms_O.get_total_energy()
                        O_ads = energy_O - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
                        df.at[dz, spin] = formation_energy

                    if os.path.exists(atoms_path):
                        atoms = read(atoms_path)
                        energy = atoms.get_total_energy()
                        # numb_N = len([atom for atom in atoms if atom.symbol == 'N'])
                        # numb_C = len([atom for atom in atoms if atom.symbol == 'C'])
                        # formation_energy = energy - metal_df.at[metal, 'energy'] - numb_C * carbon - numb_N * nitrogen
                        OH_ads = energy - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
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

                relaxed_path = os.path.join(path, 'relaxed_', 'moments.json')
                if os.path.exists(relaxed_path):
                    atoms = read(relaxed_path)
                    zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                    zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                    dz_relaxed = abs(zN - zM)
                    energy = atoms.get_total_energy()
                    # numb_N = len([atom for atom in atoms if atom.symbol == 'N'])
                    # numb_C = len([atom for atom in atoms if atom.symbol == 'C'])
                    # formation_energy = energy - metal_df.at[metal, 'energy'] - numb_C * carbon - numb_N * nitrogen
                    formation_energy = energy - metal_df.at[metal, 'energy'] - 26 * carbon - 4 * nitrogen
                    df_relaxed.at[dz_relaxed, spin] = formation_energy

                    try:
                        magmoms = atoms.get_magnetic_moments()
                        for atom in atoms:
                            if atom.symbol not in ['N', 'C', 'O', 'H']:
                                df_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                    except:
                        df_relaxed_mag.at[dz_relaxed, spin] = 0

        if 'HS' in df.columns and 'LS' in df.columns:
            df_rel['HS-LS'] = df['HS'] - df['LS']
        else:
            print(f"Warning: 'HS' or 'LS' column not found in df for {row_key}_{metal}")

        if 'HS' in df_relaxed.columns and 'LS' in df_relaxed.columns:
            df_relaxed_rel['HS-LS'] = df_relaxed['HS'] - df_relaxed['LS']
        else:
            print(f"Warning: 'HS' or 'LS' column not found in df_relaxed for {row_key}_{metal}")

        # combining(df=df, df_relaxed=df_relaxed, tsv_filename=tsv_filename)
        # combining(df=df_rel, df_relaxed=df_relaxed_rel, tsv_filename=tsv_rel_filename)
        # combining(df=df_mag, df_relaxed=df_relaxed_mag, tsv_filename=tsv_mag_filename)

        # plotting(df=df, df_relaxed=df_relaxed, dzs=dzs, spins=spins, 
        #          ylabel='Formation energy (eV)', png_filename=png_filename)
        plotting(df=df_rel, df_relaxed=df_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                 ylabel='Spin crossover energy (eV)', png_filename=png_rel_filename)
        plotting(df=df_mag, df_relaxed=df_relaxed_mag, dzs=dzs, spins=spins, ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                 ylabel='Magnetic Moments', png_filename=png_mag_filename)

def combining(df, df_relaxed, tsv_filename):
    combined_df = pd.concat([df, df_relaxed])
    combined_df.to_csv(tsv_filename, sep='\t', float_format='%.2f')
    print(f"Data saved to {tsv_filename}")

def plot_smooth_line(x, y, color, label):
    x_new = np.linspace(min(x), max(x), 300)
    spl = make_interp_spline(x, y, k=3)  # Smoothing spline
    y_smooth = spl(x_new)
    plt.plot(x_new, y_smooth, color=color, label=label)
    plt.scatter(x, y, color=color)  # Add markers without label

def plotting(df, df_relaxed, dzs, spins, ylabel, png_filename, ymin=None, ymax=None, yticks=None, color=None):
    plt.figure(figsize=(8, 6))
    for column in df.columns:
        filtered_df = df[column].dropna()
        if not filtered_df.empty:
            x = filtered_df.index
            y = filtered_df.values
            if color:
                plot_smooth_line(x, y, color, f'{column} (fixed)')
            else:
                plot_smooth_line(x, y, spins[column], f'{column} (fixed)')
    for column in df_relaxed.columns:
        filtered_df = df_relaxed[column].dropna()
        if not filtered_df.empty:
            x = filtered_df.index
            y = filtered_df.values
            if color:
                plt.scatter(x, y, marker='x', color=color, label=f'{column} (relaxed)')
            else:
                plt.scatter(x, y, marker='x', color=spins[column], label=f'{column} (relaxed)')
    if color:
        plt.axhline(y=0.0, color='blue', linestyle='--')
        plt.axhline(y=0.8, color='red', linestyle='--')
    plt.xticks(dzs)
    plt.xlabel('dz')
    plt.ylabel(ylabel)
    if ymin and ymax:
        plt.ylim(ymin, ymax)
    if yticks is not None:
        plt.yticks(yticks)
    plt.legend()
    plt.tight_layout()
    plt.savefig(png_filename, bbox_inches="tight")
    print(f"Figure saved as {png_filename}")
    plt.close()

if __name__ == '__main__':
    main()