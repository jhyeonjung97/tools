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

rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
specific_metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
colors = {'3d': 'blue', '4d': 'green', '5d': 'red'}
dirs = ['relaxed', 'oh', 'o', 'ooh']
adsorbates = ['clean', 'OH', 'O', 'OOH']

def main():
    df = pd.DataFrame()
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals, start=2):
            for d, dir in enumerate(dirs):
                ads = adsorbates[d]
                df.at[metal, 'row'] = row_key
                df.at[metal, 'color'] = colors[row_key]
                path = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/{m}_{metal}/most_stable/{dir}'
                path_DONE = os.path.join(path, 'DONE')
                path_json = os.path.join(path, 'final_with_calculator.json')
                if os.path.exists(path_DONE) and os.path.exists(path_json):
                    atoms = read(path_json)
                    energy = atoms.get_total_energy()
                    df.at[metal, ads] = energy
                else:
                    df.at[metal, ads] = np.nan
    for m, metal in enumerate(specific_metals, start=1):
        for i, ads in enumerate(['O', 'OH', 'OOH'], start=1):
            path = f'/pscratch/sd/j/jiuy97/6_MNC/{i}_{ads}/{m}_{metal}/most_stable/relaxed'
            path_DONE = os.path.join(path, 'DONE')
            path_json = os.path.join(path, 'final_with_calculator.json')
            if os.path.exists(path_DONE) and os.path.exists(path_json):
                atoms = read(path_json)
                energy = atoms.get_total_energy()
                if df.loc[metal, ads] > energy:
                    df.at[metal, ads] = energy
                
    df['G_'] = df['clean']
    df['G_OH'] = df['OH'] + OH_corr
    df['G_O'] = df['O'] + O_corr
    df['G_OOH'] = df['OOH'] + OOH_corr
    
    df['dG_OH'] = df['G_OH'] - df['G_'] - hydroxide_G
    df['dG_O'] = df['G_O'] - df['G_'] - oxygen_G 
    df['dG_OOH'] = df['G_OOH'] - df['G_'] - oxygen_G - hydroxide_G

    scaling('dG_OH', 'dG_O', 'OH', 'O', df)
            # xmin=-2.5, xmax=1.5, ymin=-4.5, ymax=4.5, txtx=0.1, txty=0.8)
    scaling('dG_OH', 'dG_OOH', 'OH', 'OOH', df)
            # xmin=0.0, xmax=1.5, ymin=3.5, ymax=4.5, txtx=0.1, txty=0.2)

    df = df.drop(columns='color')
    df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/contour/scaling_relationship_full.csv', sep=',') #, float_format='%.2f')
    df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/contour/scaling_relationship_full.tsv', sep='\t', float_format='%.2f')

def scaling(dG1, dG2, ads1, ads2, df): #, xmin, xmax, ymin, ymax):
    x, y = df[dG1], df[dG2]
    mask = np.isfinite(x) & np.isfinite(y)
    x, y, c = x[mask], y[mask], df['color'][mask]
    xx = np.linspace(min(x), max(x), 100)
    plt.figure(figsize=(4, 3), dpi=300)
    plt.scatter(x, y, c=c, s=20, zorder=4)
    for xi, yi, metal in zip(x, y, df.index):
        plt.annotate(f'{metal}', (float(xi), float(yi)), textcoords="offset points", 
                     xytext=(0, 5), ha='center', color='black', fontsize=8, zorder=5)
    coeffs = np.polyfit(x, y, 1)
    line = np.poly1d(coeffs)
    plt.plot(xx, line(xx), label=fr'$\Delta$G$_{{\sf {ads2}}}$ (trend)', linestyle='-', color='black', zorder=3)
    if ads1 == 'OH' and ads2 == 'O':
        plt.plot(xx, 2.00*xx-0.75, linestyle='--', color='lightgray', zorder=0)
        plt.plot(xx, 1.64*xx+0.95, linestyle='--', color='darkgray', zorder=1)
        plt.plot(xx, 1.49*xx+1.21, linestyle='--', color='lightblue', zorder=2)
        plt.text(0.34, 0.23, 'y = 2.00x - 0.75 (metal)', transform=plt.gca().transAxes, fontsize=10, color='lightgray')
        plt.text(0.34, 0.17, 'y = 1.64x + 0.95 (oxide)', transform=plt.gca().transAxes, fontsize=10, color='darkgray')
        plt.text(0.15, 0.05, 'y = 1.49x + 1.21 ((oxi)(hydro)oxides)', transform=plt.gca().transAxes, fontsize=10, color='lightblue')
    elif ads1 == 'OH' and ads2 == 'OOH':
        plt.plot(xx, 1.06*xx+3.16, linestyle='--', color='lightgray', zorder=0)
        plt.plot(xx, 1.05*xx+3.01, linestyle='--', color='darkgray', zorder=1)
        plt.plot(xx, 0.84*xx+3.14, linestyle='--', color='lightblue', zorder=2)
        plt.text(0.34, 0.23, 'y = 1.06x + 3.16 (metal)', transform=plt.gca().transAxes, fontsize=10, color='lightgray')
        plt.text(0.34, 0.17, 'y = 1.05x + 3.01 (oxide)', transform=plt.gca().transAxes, fontsize=10, color='darkgray')
        plt.text(0.15, 0.05, 'y = 0.84x + 3.14 ((oxi)(hydro)oxides)', transform=plt.gca().transAxes, fontsize=10, color='lightblue')
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    y_pred = line(x)
    ss_total = np.sum((y - np.mean(y))**2)
    ss_residual = np.sum((y - y_pred)**2)
    r_squared = 1 - (ss_residual / ss_total)
    equation_with_r2 = f'{equation} (R² = {r_squared:.2f})'
    plt.text(0.34, 0.11, equation_with_r2, transform=plt.gca().transAxes, fontsize=10, color='black')
    plt.xlabel(fr'$\Delta$G$_{{\sf {ads1}}}$ (eV)', fontsize='large')
    plt.ylabel(fr'$\Delta$G$_{{\sf {ads2}}}$ (eV)', fontsize='large')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='3d', markerfacecolor='blue', markersize=7),
        Line2D([0], [0], marker='o', color='w', label='4d', markerfacecolor='green', markersize=7),
        Line2D([0], [0], marker='o', color='w', label='5d', markerfacecolor='red', markersize=7)
    ]
    plt.legend(handles=legend_elements, loc='upper left')
    # plt.xlim(xmin, xmax)
    # plt.ylim(ymin, ymax)
    plt.tight_layout()
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/contour/scaling_relationship_full_{ads1}_{ads2}.png')
    print(f"Figure saved as scaling_relationship_full_{ads1}_{ads2}.png")
    plt.close()
    
if __name__ == '__main__':
    
    # Water and hydrogen properties
    water_E = -14.23091949
    water_Cv = 0.103
    water_TS = 0.675
    water_ZPE = 0.558
    water_G = water_E + water_Cv - water_TS + water_ZPE
    
    hydrogen2_E = -6.77149190
    hydrogen2_Cv = 0.0905
    hydrogen2_TS = 0.408
    hydrogen2_ZPE = 0.268
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
    
