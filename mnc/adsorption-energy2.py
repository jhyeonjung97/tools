import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np
from matplotlib.patches import Wedge

def main():
    rows = ['3d', '3d', '3d', '3d', '4d', '5d']
    numbs = [5, 6, 7, 8, 4, 4]
    metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
    adsorbates = ['clean', 'O', 'OH', 'OOH']
    ymin = [-1, -1, -1, -1, -4, -4]
    ymax = [5, 5, 5, 5, 3, 3]
    dzs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    colors = {'LS': '#ff7f0e', 'IS': '#279ff2', 'HS': '#9467bd'}
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
        
        df_energy = pd.DataFrame()
        df_spin = pd.DataFrame()
        for a, ads in enumerate(adsorbates):
            if row != '3d' and ads == 'OOH':
                df_energy[ads] = np.nan
                df_spin[ads] = np.nan
                continue
                
            energy_path = f'/pscratch/sd/j/jiuy97/6_MNC/figures/formation_energy/{row}_{numb}{metal}_{ads}.tsv'
            energy_data = pd.read_csv(energy_path, sep=',')
            ms_columns = [col for col in energy_data.columns if 'MS' in col]
            if ms_columns:
                energy_data['MS'] = energy_data[ms_columns].bfill(axis=1).iloc[:, 0]
            else:
                energy_data['MS'] = np.nan
            df_energy[ads] = energy_data['MS']
            
            if ads == 'clean':
                spin_path = f'/pscratch/sd/j/jiuy97/6_MNC/{a}_{ads}/{row}/{numb}_{metal}/lowest.tsv'
            else:
                spin_path = f'/pscratch/sd/j/jiuy97/6_MNC/{a}_{ads}/{m+1}_{metal}/lowest.tsv'
            spin_data = pd.read_csv(spin_path, sep=',')
            df_spin[ads] = spin_data['spin_state']

        df_energy = df_energy.head(7)
        df_energy.index = dzs
        df_spin.index = dzs
        
        df_energy['OHads'] = (df_energy['OH'] + OH_corr) - df_energy['clean'] - hydroxide_G
        df_energy['Oads'] = (df_energy['O'] + O_corr) - df_energy['clean'] - oxygen_G
        df_energy['OOHads'] = (df_energy['OOH'] + OOH_corr) - df_energy['clean'] - oxygen_G - hydroxide_G
        
        tsv_filename = os.path.join(save_path, f"adsorption_energy_{m+1}{metal}.tsv")
        df_energy.to_csv(tsv_filename, sep=',', float_format='%.2f')
        print(f"Data saved as adsorption_energy_{m+1}{metal}.tsv")
        
        fig, ax = plt.subplots(figsize=(3, 5), dpi=300)
        for ads in ['OH', 'O', 'OOH']:
            if row != '3d' and ads == 'OOH':
                continue
            column = f'{ads}ads'
            x = df_energy.index
            y = df_energy[column]
            plt.plot(x, y, color='black', label=f'E_{ads}',lw=0.75, zorder=1)
            for xi, yi in zip(x, y):
                left_wedge = Wedge((xi, yi), 0.05, 90, 270, lw=0.75, zorder=2, 
                                   facecolor=colors[df_spin.loc[xi, 'clean']], edgecolor='black')
                right_wedge = Wedge((xi, yi), 0.05, 270, 90, lw=0.75, zorder=2, 
                                    facecolor=colors[df_spin.loc[xi, ads]], edgecolor='black')
                ax.add_patch(left_wedge)
                ax.add_patch(right_wedge)
        plt.xlabel('Δz (Å)')
        plt.ylabel('Adsorption energy (eV)')
        # plt.xlim(0.0, 1.2)
        plt.ylim(ymin[m], ymax[m])
        plt.yticks(np.arange(ymin[m], ymax[m] + 0.1, 0.5))
        plt.xticks(dzs)
        plt.tight_layout()
        png_filename = os.path.join(save_path, f"adsorption_energy_{m+1}{metal}.png")
        plt.savefig(png_filename, bbox_inches="tight")
        print(f"Figure saved as adsorption_energy_{m+1}{metal}.png")
        plt.close()
    
if __name__ == '__main__':
    main()
