import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np

def main():
    rows = ['3d', '3d', '3d', '3d', '4d', '5d']
    numbs = [5, 6, 7, 8, 4, 4]
    metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
    adsorbates = ['clean', 'O', 'OH']
    save_path = '/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix'

    # Water and hydrogen properties
    water_E = -14.23797429
    water_Cv = 0.103
    water_TS = 0.675
    water_ZPE = 0.558
    water_G = water_E + water_Cv - water_TS + water_ZPE

    hydrogen2_E = -6.77412273
    hydrogen2_Cv = 0.0905
    hydrogen2_TS = 0.408
    hydrogen2_ZPE = 0.273
    hydrogen2_G = hydrogen2_E + hydrogen2_Cv - hydrogen2_TS + hydrogen2_ZPE
    hydrogen_G = hydrogen2_G / 2

    oxygen_G = water_G - hydrogen2_G
    hydroxide_G = water_G - hydrogen_G

    # Correction values
    OH_Cv = 0.042
    OH_TS = 0.066
    OH_ZPE = 0.376
    OH_corr = OH_Cv - OH_TS + OH_ZPE

    O_Cv = 0.034
    O_TS = 0.060
    O_ZPE = 0.064
    O_corr = O_Cv - O_TS + O_ZPE

    OOH_Cv = 0.077
    OOH_TS = 0.134
    OOH_ZPE = 0.471
    OOH_corr = OOH_Cv - OOH_TS + OOH_ZPE
    
    for m, metal in enumerate(metals):
        row = rows[m]
        numb = numbs[m]
        
        df = pd.DataFrame()
        for ads in adsorbates:
            energy_path = f'/pscratch/sd/j/jiuy97/6_MNC/figures/formation_energy/{row}_{numb}{metal}_{ads}.tsv'
            energy = pd.read_csv(energy_path, sep='\t')
            ms_columns = [col for col in energy.columns if 'MS' in col]
            if ms_columns:
                energy['MS'] = energy[ms_columns].bfill(axis=1).iloc[:, 0]
            df[ads] = energy['MS']

        df = df.head(7)
        df.index = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
        
        df['Oads'] = (df['O'] + O_corr) - df['clean'] - oxygen_G
        df['OHads'] = (df['OH'] + OH_corr) - df['clean'] - hydroxide_G
        
        tsv_filename = os.path.join(save_path, f"{m+1}{metal}_{ads}ads.tsv")
        df.to_csv(tsv_filename, sep='\t', float_format='%.2f')
        print(f"Data saved as {m+1}{metal}_Eads.tsv")

        plt.figure(figsize=(4, 3), dpi=300)
        plt.plot(df.index, df['Oads'], marker='o', label='E_O')
        plt.plot(df.index, df['OHads'], marker='o', label='E_OH')

        plt.xlabel('dz (A)')
        plt.ylabel('Adsorption energy (eV)')
        plt.legend()
        plt.tight_layout()
        png_filename = os.path.join(save_path, f"{m+1}{metal}_{ads}ads.png")
        plt.savefig(png_filename, bbox_inches="tight")
        print(f"Figure saved as {m+1}{metal}_Eads.png")
        plt.close()        
        
if __name__ == '__main__':
    main()
