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
filename = 'heo_relative_energy'
ref_filename = 'heo_references'
chg_filename = 'heo_bader_charge'
mag_filename = 'heo_magnetic_moments'
occ_filename = 'heo_eg_occupancies'
gap_filename = 'heo_band_gap'
dos_filename = 'heo_density_of_states'

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
                            
    saving(df, filename)
    saving(df_chg, chg_filename)
    saving(df_mag, mag_filename)
    saving(df_occ, occ_filename)
    saving(df_ref, ref_filename)

    for i in range(5):
        df_ref.at[i, 'energy'] = 0
        
    plotting(pattern='energy', xlabel='Relative energy (eV)', filename=filename, figsize=(6, 6), min=-0.4, max=+0.4, tick=0.08) 
    plotting(pattern='bandgap', xlabel='Band gap (eV)', filename=gap_filename, figsize=(10, 6), min=+0.0, max=+2.8, tick=0.20)
    plotting(pattern='Md2Op', xlabel='M3d - O2p (eV)', filename=dos_filename, figsize=(8, 6), min=+0.0, max=+2.0, tick=0.20)
    
    plotting_adv(df=df_mag, df_ref=df_ref, pattern='magmom', xlabel='Magnetic moments (uB)', filename=mag_filename, min=0, max=5, tick=0.1)
    plotting_adv(df=df_chg, df_ref=df_ref, pattern='charge', xlabel='Bader charge (e-)', filename=chg_filename, min=1.0, max=2.1, tick=0.02)
    plotting_adv(df=df_occ, df_ref=df_ref, pattern='eg_occ', xlabel='e_g occupancy (e-)', filename=occ_filename, min=2.2, max=4.0, tick=0.04)
    
def saving(df, filename):
    df.to_csv(f'{filename}.tsv', sep='\t', float_format='%.2f')
    print(f"Data saved to {filename}.tsv")

def plotting(pattern, xlabel, filename, min, max, tick):
    
    bins=np.arange(min, max, tick*1); width=tick*0.9; xticks=np.arange(min, max, tick*2); xmin=min-tick*5; xmax=max+tick*5
    
    plt.figure(figsize=figsize)
    plt.hist(df[pattern].dropna(), bins=bins, alpha=0.5, width=width)
    for i in range(5):
        plt.axvline(x=df_ref.at[i, pattern], color=clrs[i], linestyle='--')
    plt.axvline(x=0, color='gray', linestyle='--')
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.xticks(xticks)
    plt.xlim(xmin, xmax)
    plt.savefig(f'{filename}.png', bbox_inches="tight")
    print(f"Figure saved as {filename}")
    plt.close()

def plotting_adv(df, df_ref, pattern, xlabel, filename, min, max, tick,
                 figsize1=(8, 6), figsize2=(12, 6)):

    bins1=np.arange(min, max, tick*1); width1=tick*0.9; xticks1=np.arange(min, max, tick*15); xmin1=min-tick*5; xmax1=max+tick*5
    bins2=np.arange(min, max, tick*2); width2=tick*2.0; xticks2=np.arange(min, max, tick*15); xmin2=min-tick*5; xmax2=max+tick*5

    for i, column in enumerate(df_chg.columns):
        plt.figure(figsize=figsize1)
        plt.hist(df_mag[column].dropna(), bins=bins1, alpha=0.5, color=clrs[i], label=str(column), width=width1)
        plt.axvline(x=df_ref.at[i, pattern], color=clrs[i], linestyle='--')
        plt.xlabel(xlabel)
        plt.ylabel('Frequency')
        plt.xticks(xticks1)
        plt.xlim(xmin1, xmax1)
        plt.legend(title="B sites")
        plt.savefig(f'{filename}_{column}.png', bbox_inches="tight")
        print(f"Figure saved as {filename}_{column}.png")
        plt.close()

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