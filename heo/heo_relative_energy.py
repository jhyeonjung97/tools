import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import pandas as pd
import numpy as np
import os
import re

# Define the metals and initialize dataframes
prvs = ['Cr', 'Mn', 'Fe', 'Co', 'Ni']
clrs = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
df = pd.DataFrame()
df_chg = pd.DataFrame()
df_mag = pd.DataFrame()
df_ref = pd.DataFrame()
df_occ = pd.DataFrame()
numb = [0] * 5

# Filenames for saving the data and plots
tsv_filename = 'heo_relative_energy.tsv'
png_filename = 'heo_relative_energy.png'


ref_filename = 'heo_references'
chg_filename = 'heo_bader_charge'
mag_filename = 'heo_magnetic_moments'
occ_filename = 'heo_eg_occupancies'

png_gap_filename = 'heo_band_gap.png'
png_dos_filename = 'heo_density_of_states.png'

# Regular expressions to extract band gap and DOS information
pattern_gap = re.compile(r"Band Gap:\s+([\d.]+)\s+eV")
pattern_dos = re.compile(r"Average Energy \(band center\):\s+([-+]?\d*\.\d+|\d+)")

def main():
    for i in range(5):
        path = f'/scratch/x2755a09/4_HEO/pure/{i+1}_{prvs[i]}/moments.json'
        chg_path = f'/scratch/x2755a09/4_HEO/pure/{i+1}_{prvs[i]}/atoms_bader_charge.json'
        gap_path = f'/scratch/x2755a09/4_HEO/pure/{i+1}_{prvs[i]}/gap.txt'
        dos_path = f'/scratch/x2755a09/4_HEO/pure/{i+1}_{prvs[i]}/dos.txt'
        occ_path = f'/scratch/x2755a09/4_HEO/pure/{i+1}_{prvs[i]}/occ.tsv'
        
        if os.path.exists(path):
            atoms = read(path)
            magmoms = atoms.get_magnetic_moments()
            df_ref.at[i, 'energy'] = atoms.get_total_energy()
            df_ref.at[i, 'magmom'] = mean([abs(magmoms[atom.index]) for atom in atoms if atom.symbol == prvs[i]])
        if os.path.exists(chg_path):
            atoms = read(chg_path)
            charges = atoms.get_initial_charges()
            df_ref.at[i, 'charge'] = mean([abs(charges[atom.index]) for atom in atoms if atom.symbol == prvs[i]])
        if os.path.exists(gap_path):
            with open(gap_path, 'r') as file:
                lines = file.read()
                match = pattern_gap.search(lines)
                if match:
                    df_ref.at[i, 'bandgap'] = float(match.group(1))
        if os.path.exists(dos_path):
            with open(dos_path, 'r') as file:
                lines = file.read()
                matches = pattern_dos.findall(lines)
                if len(matches) == 2:
                    df_ref.at[i, 'Md2Op'] = float(matches[0]) - float(matches[1])
        if os.path.exists(occ_path):
            df_occ_tmp = pd.read_csv(occ_path, delimiter='\t', index_col=0)
            df_ref.at[i, 'eg_occ'] = df_occ_tmp[['occ4', 'occ5', 'occ9', 'occ10']].sum(axis=1).mean()
    
    for i in range(60):
        path = f'/scratch/x2755a09/4_HEO/{i:02d}_/final_with_calculator.json'
        chg_path = f'/scratch/x2755a09/4_HEO/{i:02d}_/atoms_bader_charge.json'
        gap_path = f'/scratch/x2755a09/4_HEO/{i:02d}_/gap.txt'
        dos_path = f'/scratch/x2755a09/4_HEO/{i:02d}_/dos.txt'
        occ_path = f'/scratch/x2755a09/4_HEO/{i:02d}_/occ.tsv'

        indice = {metal: [] for metal in prvs}
        if os.path.exists(path):
            atoms = read(path)
            energy = atoms.get_total_energy()
            magmoms = atoms.get_magnetic_moments()
            for m, metal in enumerate(prvs):
                numb[m] = len([atom for atom in atoms if atom.symbol == metal])
                for atom in atoms:
                    if atom.symbol == metal:
                        indice[metal].append(atom.index)
                df_mag.at[i, metal] = mean([abs(magmoms[idx]) for idx in indice[metal]])
            relative_energy = energy - sum(numb[m] * df_ref.at[m, 'energy'] / 8 for m, metal in enumerate(prvs))
            df.at[i, 'energy'] = relative_energy
        if os.path.exists(path):
            atoms = read(chg_path)
            charges = atoms.get_initial_charges()
            for metal in prvs:
                df_chg.at[i, metal] = mean([abs(charges[atom.index]) for atom in atoms if atom.symbol == metal])
        if os.path.exists(gap_path):
            with open(gap_path, 'r') as file:
                lines = file.read()
                match = pattern_gap.search(lines)
                if match:
                    df.at[i, 'bandgap'] = float(match.group(1))
        if os.path.exists(dos_path):
            with open(dos_path, 'r') as file:
                lines = file.read()
                matches = pattern_dos.findall(lines)
                if len(matches) == 2:
                    df.at[i, 'Md2Op'] = float(matches[0]) - float(matches[1])
        if os.path.exists(occ_path):                
            df_occ_tmp = pd.read_csv(occ_path, delimiter='\t', index_col=0)
            for metal in prvs:
                tmp = []
                for idx in indice[metal]:
                    for o in [4, 5, 9, 10]:
                        tmp.append(df_occ_tmp.loc[f'atom_{idx+1}', f'occ{o}'])
                if tmp:
                    df_occ.at[i, metal] = mean(tmp)*4
                else:
                    df_occ.at[i, metal] = np.nan  # Handle case where tmp is empty
                            
    saving(df, tsv_filename)
    saving(df_chg, tsv_chg_filename)
    saving(df_mag, tsv_mag_filename)
    saving(df_occ, tsv_occ_filename)
    saving(df_ref, tsv_ref_filename)

    for i in range(5):
        df_ref.at[i, 'energy'] = 0
        
    plotting(pattern='energy', xlabel='Relative energy (eV)', filename=png_filename, 
             figsize=(6, 6), bins=np.arange(-2.0, 0.0, 0.1), width=0.09, xticks=np.arange(-2.0, 0.1, 0.2), xmin=-1.5, xmax=0.1)
    plotting(pattern='bandgap', xlabel='Band gap (eV)', filename=png_gap_filename, 
             figsize=(10, 6), bins=np.arange(0.0, 2.2, 0.1), width=0.09, xticks=np.arange(0.0, 2.9, 0.2), xmin=-0.1, xmax=2.9)
    plotting(pattern='Md2Op', xlabel='M3d - O2p (eV)', filename=png_dos_filename, 
             figsize=(8, 6), bins=np.arange(0.4, 2.8, 0.1), width=0.09, xticks=np.arange(0.0, 2.3, 0.2), xmin=-0.1, xmax=2.3)
    
    plotting_adv(df=df_mag, df_ref=df_ref, pattern='magmom', xlabel='Magnetic moments (uB)', filename=mag_filename,
                 bins1=np.arange(0, 6, 0.1), width1=0.09, xticks1=np.arange(0, 6, 1), xmin1=-0.5, xmax1=5.5, 
                 bins2=np.arange(0, 6, 0.2), width2=0.2, xticks2=np.arange(0, 6, 1), xmin2=-0.5, xmax2=5.5)
    plotting_adv(df=df_chg, df_ref=df_ref, pattern='charge', xlabel='Bader charge (e-)', filename=chg_filename,
                 bins1=np.arange(1.0, 2.1, 0.02), width1=0.018, xticks1=np.arange(1.0, 2.1, 0.1), xmin1=0.95, xmax1=2.05, 
                 bins2=np.arange(1.0, 2.1, 0.04), width2=0.04, xticks2=np.arange(1.0, 2.1, 0.1), xmin2=0.95, xmax2=2.05)
    plotting_adv(df=df_occ, df_ref=df_ref, pattern='eg_occ', xlabel='e_g occupancy (e-)', filename=occ_filename,
                 bins1=np.arange(2.2, 4.0, 0.04), width1=0.036, xticks1=np.arange(2.2, 4.2, 0.2), xmin1=2.1, xmax1=4.1, 
                 bins2=np.arange(2.2, 4.0, 0.08), width2=0.08, xticks2=np.arange(2.2, 4.2, 0.2), xmin2=2.1, xmax2=4.1)
    
def saving(df, filename):
    df.to_csv(f'{filename}.tsv', sep='\t', float_format='%.2f')
    print(f"Data saved to {filename}.tsv")

def plotting(pattern, xlabel, filename, 
             figsize, bins, width, xticks, xmin, xmax):
    plt.figure(figsize=figsize)
    plt.hist(df[pattern].dropna(), bins=bins, alpha=0.5, width=width)
    for i in range(5):
        plt.axvline(x=df_ref.at[i, pattern], color=clrs[i], linestyle='--')
    plt.axvline(x=0, color='gray', linestyle='--')
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.xticks(xticks)
    plt.xlim(xmin, xmax)
    plt.savefig(filename, bbox_inches="tight")
    print(f"Figure saved as {filename}")
    plt.close()

def plotting_adv(df, df_ref, pattern, xlabel, filename,
                 bins1, width1, xticks1, xmin1, xmax1,
                 bins2, width2, xticks2, xmin2, xmax2, 
                 figsize1=(8, 6), figsize2=(12, 6)):
    # for i, column in enumerate(df_chg.columns):
    #     plt.figure(figsize=figsize1)
    #     plt.hist(df_mag[column].dropna(), bins=bins1, alpha=0.5, color=clrs[i], label=str(column), width=width1)
    #     plt.axvline(x=df_ref.at[i, pattern], color=clrs[i], linestyle='--')
    #     plt.xlabel(xlabel)
    #     plt.ylabel('Frequency')
    #     plt.xticks(xticks1)
    #     plt.xlim(xmin1, xmax1)
    #     plt.legend(title="B sites")
    #     plt.savefig(f'{filename}_{column}.png', bbox_inches="tight")
    #     print(f"Figure saved as {filename}_{column}.png")
    #     plt.close()

    plt.figure(figsize=figsize2)
    for i in range(5):
        plt.axvline(x=df_ref.at[i, pattern], color=clrs[i], linestyle='--')
    bins = bins2
    bin_width = width2 / (len(df.columns) + 1)  # Calculate new width for each bar
    for idx, column in enumerate(df.columns):
        plt.hist(df[column].dropna(), bins=bins2 + idx * bin_width, alpha=0.5, label=str(column), width=bin_width)
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.xticks(xticks2)
    plt.xlim(xmin2, xmax2)
    plt.legend(title="B sites")
    plt.savefig(f'{filename}.png', bbox_inches="tight")
    print(f"Figure saved as {filename}.png")
    plt.close()

if __name__ == '__main__':
    main()