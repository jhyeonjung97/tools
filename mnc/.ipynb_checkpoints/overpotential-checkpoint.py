import os
import glob
import numpy as np
import pandas as pd
from ase.io import read
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import make_interp_spline

OERs = ['H2O->*OH','*OH->*O', '*O->*OOH', '*OOH->O2']
ORRs = ['O2->*OOH','*OOH->*O', '*O->*OH', '*OH->H2O']
rows = ['3d', '3d', '3d', '3d', '4d', '5d']
groups = ['5', '6', '7', '8', '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
adsorbates = ['clean', 'O', 'OH']
spins = ['LS', 'IS', 'HS']
colors = ['#ff7f0e', '#279ff2', '#9467bd']

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

root = '/pscratch/sd/j/jiuy97/6_MNC/figure'

def main():
    for m, metal in enumerate(metals):
        row = rows[m]
        group = groups[m]
        energies = {}
        gibbs_energies = pd.DataFrame()
        spin_cross_over = pd.DataFrame()
        for adsorbate in adsorbates:
            tsv_path = os.path.join(root, f'{row}_{group}{metal}_{adsorbate}.tsv')
            energies[adsorbate] = pd.read_csv(tsv_path, sep='\t', index_col=0)
            energies[adsorbate] = energies[adsorbate].head(7)
            # for spin in spins:
            #     if spin in energies[adsorbate].columns:
            #         energies[adsorbate] = energies[adsorbate].drop(columns=[spin])
            energies[adsorbate]['energy'] = energies[adsorbate].min(axis=1, skipna=True)
            energies[adsorbate]['spin'] = energies[adsorbate].apply(
                lambda row: row.idxmin(skipna=True) if row.notna().any() else None, axis=1)
            energies[adsorbate]['spin'] = energies[adsorbate]['spin'].apply(
                lambda x: f'MS({x})' if x in ['LS', 'IS', 'HS'] else x)
            
        gibbs_energies['G_'] = energies['clean']['energy']
        gibbs_energies['G_OH'] = energies['OH']['energy'] + OH_corr
        gibbs_energies['G_O'] = energies['O']['energy'] + O_corr
        
        G_ = min(gibbs_energies['G_'])
        G_OH = min(gibbs_energies['G_OH'])
        G_O= min(gibbs_energies['G_O'])
        
        gibbs_energies['dG_OH'] = gibbs_energies['G_OH'] - gibbs_energies['G_'] - hydroxide_G
        gibbs_energies['dG_O'] = gibbs_energies['G_O'] - gibbs_energies['G_'] - oxygen_G
        gibbs_energies['dG_OOH'] = gibbs_energies['dG_OH'] + 3.2
        
        dG_OH = G_OH - G_ - hydroxide_G
        dG_O = G_O - G_ - oxygen_G
        dG_OOH = dG_OH + 3.2
        
        gibbs_energies['dG1'] = gibbs_energies['dG_OH']
        gibbs_energies['dG2'] = gibbs_energies['dG_O'] - gibbs_energies['dG_OH']
        gibbs_energies['dG3'] = gibbs_energies['dG_OOH'] - gibbs_energies['dG_O']
        gibbs_energies['dG4'] = 4.92 - gibbs_energies['dG_OOH']
        
        dG1 = dG_OH
        dG2 = dG_O - dG_OH
        dG3 = dG_OOH - dG_O
        dG4 = 4.92 - dG_OOH
        
        if gibbs_energies[['dG1', 'dG2', 'dG3', 'dG4']].notna().all().all():
            gibbs_energies['OER'] = gibbs_energies[['dG1', 'dG2', 'dG3', 'dG4']].max(axis=1) - 1.23
            gibbs_energies['ORR'] = 1.23 - gibbs_energies[['dG1', 'dG2', 'dG3', 'dG4']].min(axis=1)
        else:
            gibbs_energies['OER'] = None
            gibbs_energies['ORR'] = None
            
        if all([dG1 is not None, dG2 is not None, dG3 is not None, dG4 is not None]):
            OER = max(dG1, dG2, dG3, dG4) - 1.23
            ORR = 1.23 - min(dG1, dG2, dG3, dG4)
        else:
            OER = None
            ORR = None
            
        # gibbs_energies = gibbs_energies.round(2)
        gibbs_energies = gibbs_energies.set_index(energies['clean'].index)
            
        for index in energies['clean'].index:
            spin_cross_over.loc[index, OERs[0]] = f"{energies['clean']['spin'].loc[index]}->{energies['OH']['spin'].loc[index]}"
            spin_cross_over.loc[index, OERs[1]] = f"{energies['OH']['spin'].loc[index]}->{energies['O']['spin'].loc[index]}"
            spin_cross_over.loc[index, ORRs[2]] = f"{energies['O']['spin'].loc[index]}->{energies['OH']['spin'].loc[index]}"
            spin_cross_over.loc[index, ORRs[3]] = f"{energies['OH']['spin'].loc[index]}->{energies['clean']['spin'].loc[index]}"
    
        
        gibbs_energies.to_csv(f'{row}_{group}{metal}_gibbs.tsv', sep='\t', float_format='%.2f')
        spin_cross_over.to_csv(f'{row}_{group}{metal}_spin.tsv', sep='\t')
        print(f"Data saved to {row}_{group}{metal}_gibbs.tsv and {row}_{group}{metal}_spin.tsv")
        
        plotting(gibbs_energies=gibbs_energies, spin_cross_over=spin_cross_over, row=row, group=group, metal=metal,
                 rxn='OER', overpotential=OER, ylabel='Energy (eV)')
        plotting(gibbs_energies=gibbs_energies, spin_cross_over=spin_cross_over, row=row, group=group, metal=metal,
                 rxn='ORR', overpotential=ORR, ylabel='Energy (eV)')
        
def plot_smooth_line(x, y, color):
    try:
        x_new = np.linspace(min(x), max(x), 300)
        if len(x) > 3:
            spl = make_interp_spline(x, y, k=3)  # Smoothing spline
        else:
            spl = make_interp_spline(x, y, k=2)
        y_smooth = spl(x_new)
        plt.plot(x_new, y_smooth, color=color, zorder=1)
        plt.scatter(x, y, marker='s', edgecolors=color, facecolors='white', zorder=2)
        return y_smooth
    except ValueError as e:
        print(f"Error while creating spline: {e}")
        return None

def plotting(gibbs_energies, spin_cross_over, row, group, metal, 
             rxn, overpotential, ylabel):
    if gibbs_energies.isna().all().all():
        print("dataframe contains only NaN values, skipping plot.")
        return
    plt.figure(figsize=(4, 3))
    png_filename=f'{row}_{group}{metal}_{rxn}.png'
    filtered_gibbs_energies = gibbs_energies[rxn].dropna()
    if not filtered_gibbs_energies.empty:
        x = filtered_gibbs_energies.index
        y = filtered_gibbs_energies.values
        # plt.scatter(x, y)
        try:
            x_new = np.linspace(min(x), max(x), 300)
            if len(x) > 3:
                spl = make_interp_spline(x, y, k=3)  # Smoothing spline
            else:
                spl = make_interp_spline(x, y, k=2)
            y_smooth = spl(x_new)
            plt.plot(x_new, y_smooth, zorder=1)
            plt.scatter(x, y, marker='s', zorder=2)
        except ValueError as e:
            print(f"Error while creating spline")        
    if overpotential:
        plt.axhline(y=overpotential, color='black', linestyle='--', linewidth=0.5, zorder=0)
    plt.xlabel('dz (Å)')
    plt.ylabel(ylabel)
    plt.ylim(0.4, 1.4)
    plt.yticks(np.arange(0.4, 1.4, 0.1))
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