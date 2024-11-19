import os
import re
from statistics import mean
import glob
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
from ase.io import read, write
from matplotlib import rc

# Define the rows and spins
rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
dirs = ['relaxed', 'oh', 'o', 'ooh']
adsorbates = ['clean', 'OH', 'O', 'OOH']
df = pd.DataFrame()

def main():
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals):
            for d, dir in enumerate(dirs):
                ads = adsorbates[d]
                path = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/{m+2}_{metal}/most_stable/{dir}'
                path_DONE = os.path.join(path, 'DONE')
                path_json = os.path.join(path, 'final_with_calculator.json')
                if os.path.exists(path_DONE) and os.path.exists(path_json):
                    atoms = read(path_json)
                    energy = atoms.get_total_energy()
                    df.at[metal, ads] = energy
                else:
                    df.at[metal, ads] = np.nan

    df['G_'] = df['clean']
    df['G_OH'] = df['OH'] + OH_corr
    df['G_O'] = df['O'] + O_corr
    df['G_OOH'] = df['OOH'] + OOH_corr
    
    df['dG_OH'] = df['G_OH'] - df['G_'] - hydroxide_G
    df['dG_O'] = df['G_O'] - df['G_'] - oxygen_G 
    df['dG_OOH'] = df['G_OOH'] - df['G_'] - oxygen_G - hydroxide_G

    df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/scaling_relationship_full.tsv', sep='\t', float_format='%.2f')

    
    
def scaling(x, y, ads1, ads2, scaling_relationship, metals, xmin, xmax, ymin, ymax, txtx, txty):
    xx = np.linspace(min(x), max(x), 100)
    plt.figure(figsize=(4, 3), dpi=300)
    plt.scatter(x, y, c=colors[:len(x)], s=20)
    for xi, yi, metal in zip(x, y, metals):
        plt.annotate(f'{metal}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='black')
    coeffs = np.polyfit(x, y, 1)
    line = np.poly1d(coeffs)
    plt.plot(xx, line(xx), label=fr'$\Delta$G$_{{\sf {ads2}}}$ (trend)', linestyle='-', color='black')
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    plt.text(txtx, txty if coeffs[0] > 0 else 0.1, equation, transform=plt.gca().transAxes, fontsize=10, color='black')
    plt.xlabel(fr'$\Delta$G$_{{\sf {ads1}}}$ (eV)', fontsize='large')
    plt.ylabel(fr'$\Delta$G$_{{\sf {ads2}}}$ (eV)', fontsize='large')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tight_layout()
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/scaling_relationship_{ads1}_{ads2}.png')
    print(f"Figure saved as scaling_relationship_{ads1}_{ads2}.png")
    plt.close()
    
if __name__ == '__main__':
    
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
    
    main()
    
