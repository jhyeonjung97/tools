import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np

# Define the rows and spins
rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
spins = {'LS': '#ff7f0e', 'IS': '#279ff2', 'HS': '#9467bd'}
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

metal_path = '/scratch/x2755a09/3_MNC/gas/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)

def main():
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals):
            if metal not in ['Ti', 'V', 'Cr', 'Mn']:
                df = pd.DataFrame()
                df_rel = pd.DataFrame()
                df_mag = pd.DataFrame()
                df_relaxed = pd.DataFrame()
                df_relaxed_rel = pd.DataFrame()
                df_relaxed_mag = pd.DataFrame()

                df_O = pd.DataFrame()
                df_O_rel = pd.DataFrame()
                df_O_mag = pd.DataFrame()
                df_O_relaxed = pd.DataFrame()
                df_O_relaxed_rel = pd.DataFrame()
                df_O_relaxed_mag = pd.DataFrame()

                df_OH = pd.DataFrame()
                df_OH_rel = pd.DataFrame()
                df_OH_mag = pd.DataFrame()
                df_OH_relaxed = pd.DataFrame()
                df_OH_relaxed_rel = pd.DataFrame()
                df_OH_relaxed_mag = pd.DataFrame()

                tsv_filename = f'{row_key}_{m+2}{metal}_Ef.tsv'
                png_filename = f'{row_key}_{m+2}{metal}_Ef.png'
                tsv_rel_filename = f'{row_key}_{m+2}{metal}_rel.tsv'
                png_rel_filename = f'{row_key}_{m+2}{metal}_rel.png'
                tsv_mag_filename = f'{row_key}_{m+2}{metal}_mag.tsv'
                png_mag_filename = f'{row_key}_{m+2}{metal}_mag.png'

                tsv_O_filename = f'{row_key}_{m+2}{metal}_Oads.tsv'
                png_O_filename = f'{row_key}_{m+2}{metal}_Oads.png'
                tsv_O_rel_filename = f'{row_key}_{m+2}{metal}_rel_O.tsv'
                png_O_rel_filename = f'{row_key}_{m+2}{metal}_rel_O.png'
                tsv_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.tsv'
                png_O_mag_filename = f'{row_key}_{m+2}{metal}_mag_O.png'

                tsv_OH_filename = f'{row_key}_{m+2}{metal}_OHads.tsv'
                png_OH_filename = f'{row_key}_{m+2}{metal}_OHads.png'
                tsv_OH_rel_filename = f'{row_key}_{m+2}{metal}_rel_OH.tsv'
                png_OH_rel_filename = f'{row_key}_{m+2}{metal}_rel_OH.png'
                tsv_OH_mag_filename = f'{row_key}_{m+2}{metal}_mag_OH.tsv'
                png_OH_mag_filename = f'{row_key}_{m+2}{metal}_mag_OH.png'

                for spin in spins.keys():
                    path_pattern = f'/scratch/x2755a09/3_MNC/{row_key}/*_{metal}/*_{spin}'
                    matching_paths = glob.glob(path_pattern)

                    for path in matching_paths:
                        print(path)

                        for i, dz in enumerate(dzs):
                            atoms_path = os.path.join(path, f'{i}_', 'moments.json')
                            atoms_O_path = os.path.join(path, f'{i}_/1_O', 'moments.json')
                            atoms_OH_path = os.path.join(path, f'{i}_/2_OH', 'moments.json')

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

                            if os.path.exists(atoms_O_path):
                                atoms = read(atoms_O_path)
                                energy_O = atoms.get_total_energy()
                                adsorption_energy = energy_O - energy - E_O
                                df_O.at[dz, spin] = adsorption_energy
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

                            if os.path.exists(atoms_OH_path):
                                atoms = read(atoms_OH_path)
                                energy_OH = atoms.get_total_energy()
                                adsorption_energy = energy_OH - energy - E_OH
                                df_OH.at[dz, spin] = adsorption_energy
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

                        relaxed_path = os.path.join(path, 'relaxed_', 'moments.json')
                        relaxed_O_path = os.path.join(path, 'relaxed_/1_O', 'moments.json')
                        relaxed_OH_path = os.path.join(path, 'relaxed_/2_OH', 'moments.json')

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

                        if os.path.exists(relaxed_O_path):
                            atoms = read(relaxed_O_path)
                            zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                            zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                            dz_relaxed = abs(zN - zM)
                            energy_O = atoms.get_total_energy()
                            adsorption_energy = energy_O - energy - E_O
                            df_O_relaxed.at[dz_relaxed, spin] = adsorption_energy
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_O_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                            except:
                                df_O_relaxed_mag.at[dz_relaxed, spin] = 0

                        if os.path.exists(relaxed_OH_path):
                            atoms = read(relaxed_OH_path)
                            zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
                            zM = mean([atom.z for atom in atoms if atom.symbol not in ['N', 'C', 'O', 'H']])
                            dz_relaxed = abs(zN - zM)
                            energy_OH = atoms.get_total_energy()
                            adsorption_energy = energy_OH - energy - E_OH
                            df_OH_relaxed.at[dz_relaxed, spin] = adsorption_energy
                            try:
                                magmoms = atoms.get_magnetic_moments()
                                for atom in atoms:
                                    if atom.symbol not in ['N', 'C', 'O', 'H']:
                                        df_OH_relaxed_mag.at[dz_relaxed, spin] = magmoms[atom.index]
                            except:
                                df_OH_relaxed_mag.at[dz_relaxed, spin] = 0

                relative(df, df_rel)
                relative(df_O, df_O_rel)
                relative(df_OH, df_OH_rel)
                relative(df_relaxed, df_relaxed_rel)
                relative(df_O_relaxed, df_O_relaxed_rel)
                relative(df_OH_relaxed, df_OH_relaxed_rel)

                combining(df=df, df_relaxed=df_relaxed, tsv_filename=tsv_filename)
                combining(df=df_rel, df_relaxed=df_relaxed_rel, tsv_filename=tsv_rel_filename)
                combining(df=df_mag, df_relaxed=df_relaxed_mag, tsv_filename=tsv_mag_filename)
                combining(df=df_O, df_relaxed=df_O_relaxed, tsv_filename=tsv_O_filename)
                combining(df=df_O_rel, df_relaxed=df_O_relaxed_rel, tsv_filename=tsv_O_rel_filename)
                combining(df=df_O_mag, df_relaxed=df_O_relaxed_mag, tsv_filename=tsv_O_mag_filename)
                combining(df=df_OH, df_relaxed=df_OH_relaxed, tsv_filename=tsv_OH_filename)
                combining(df=df_OH_rel, df_relaxed=df_OH_relaxed_rel, tsv_filename=tsv_OH_rel_filename)
                combining(df=df_OH_mag, df_relaxed=df_OH_relaxed_mag, tsv_filename=tsv_OH_mag_filename)

                plotting(df=df, df_relaxed=df_relaxed, dzs=dzs, spins=spins, 
                         ylabel='Formation energy (eV)', png_filename=png_filename)
                plotting(df=df_rel, df_relaxed=df_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                         ylabel='Spin crossover energy (eV)', png_filename=png_rel_filename)
                plotting(df=df_mag, df_relaxed=df_relaxed_mag, dzs=dzs, spins=spins, 
                         ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                         ylabel='Magnetic Moments', png_filename=png_mag_filename)
                plotting(df=df_O, df_relaxed=df_O_relaxed, dzs=dzs, spins=spins, 
                         ylabel='Formation energy (eV)', png_filename=png_O_filename)
                plotting(df=df_O_rel, df_relaxed=df_O_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                         ylabel='Spin crossover energy (eV)', png_filename=png_O_rel_filename)
                plotting(df=df_O_mag, df_relaxed=df_O_relaxed_mag, dzs=dzs, spins=spins, 
                         ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                         ylabel='Magnetic Moments', png_filename=png_O_mag_filename)
                plotting(df=df_OH, df_relaxed=df_OH_relaxed, dzs=dzs, spins=spins, 
                         ylabel='Formation energy (eV)', png_filename=png_OH_filename)
                plotting(df=df_OH_rel, df_relaxed=df_OH_relaxed_rel, dzs=dzs, spins=spins, color='black', 
                         ylabel='Spin crossover energy (eV)', png_filename=png_OH_rel_filename)
                plotting(df=df_OH_mag, df_relaxed=df_OH_relaxed_mag, dzs=dzs, spins=spins, 
                         ymin=-0.5, ymax=5.5, yticks=np.arange(6),
                         ylabel='Magnetic Moments', png_filename=png_OH_mag_filename)
            
def relative(df, df_rel):
    if 'HS' in df.columns and 'LS' in df.columns:
        df_rel['HS-LS'] = df['HS'] - df['LS']

def combining(df, df_relaxed, tsv_filename):
    combined_df = pd.concat([df, df_relaxed])
    combined_df.to_csv(tsv_filename, sep='\t', float_format='%.2f')
    print(f"Data saved to {tsv_filename}")

def plot_smooth_line(x, y, color, label):
    try:
        x_new = np.linspace(min(x), max(x), 300)
        if len(x) > 3:
            spl = make_interp_spline(x, y, k=3)  # Smoothing spline
        else:
            spl = make_interp_spline(x, y, k=2)
        y_smooth = spl(x_new)
        plt.plot(x_new, y_smooth, color=color, label=label)
        plt.scatter(x, y, color=color)  # Add markers without label
    except ValueError as e:
        print(f"Error while creating spline: {e}")

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