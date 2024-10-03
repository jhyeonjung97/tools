import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import make_interp_spline
from matplotlib.patches import Wedge, Rectangle

# Figure and font settings
fig_width_pt = 1.8 * 246.0
inches_per_pt = 1.0 / 72.27
golden_mean = (np.sqrt(5) - 1.0) / 2.0
fig_width = fig_width_pt * inches_per_pt
fig_height = fig_width * golden_mean
print(fig_width, fig_height)

# Constants for Gibbs energy calculations
OERs = ['H2O->*OH', '*OH->*O', '*O->*OOH', '*OOH->O2']
ORRs = ['O2->*OOH', '*OOH->*O', '*O->*OH', '*OH->H2O']
rows = ['3d', '3d', '3d', '3d', '4d', '5d']
groups = ['5', '6', '7', '8', '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
adsorbates = ['clean', 'O', 'OH']
# colors = ['#FFC3BD', '#A8E6A1', '#FFD92F', '#A0C8F8']
colors = ['blue', 'green', 'orange', 'red', 'purple', 'grey']
ms_colors = {'MS(LS)': '#ff7f0e', 'MS(IS)': '#279ff2', 'MS(HS)': '#9467bd'}

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

# Set the root directory for data files
root = '/pscratch/sd/j/jiuy97/6_MNC/figure'
relaxed_energies = {}
scaling_relationship = pd.DataFrame()

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
            relaxed_energies[adsorbate] = energies[adsorbate].iloc[7:].copy()

            for column in relaxed_energies[adsorbate].columns:
                if relaxed_energies[adsorbate][column].isna().all() and not energies[adsorbate][column].isna().all():
                    relaxed_min = energies[adsorbate][column].dropna().min()
                    relaxed_dz = energies[adsorbate][column].dropna().idxmin()
                    relaxed_energies[adsorbate].at[relaxed_dz, column] = relaxed_min
                    
            ms_columns = [col for col in relaxed_energies[adsorbate].columns if 'MS' in col]
            non_ms_columns = [col for col in relaxed_energies[adsorbate].columns if 'MS' not in col]

            if ms_columns:
                ms_data = relaxed_energies[adsorbate][ms_columns]
                if not ms_data.dropna(how='all').empty:
                    scaling_min = ms_data.min().min()
                    scaling_dz = ms_data.stack().idxmin()[0]
                else:
                    non_ms_data = relaxed_energies[adsorbate][non_ms_columns]
                    scaling_min = non_ms_data.min().min()
                    scaling_dz = non_ms_data.stack().idxmin()[0]
            else:
                non_ms_data = relaxed_energies[adsorbate][non_ms_columns]
                scaling_min = non_ms_data.min().min()
                scaling_dz = non_ms_data.stack().idxmin()[0]

            tag = '' if adsorbate == 'clean' else adsorbate
            scaling_relationship.at[metal, f'G_{tag}'] = scaling_min
            scaling_relationship.at[metal, f'dz_{tag}'] = scaling_dz

            energies[adsorbate] = energies[adsorbate].head(7)
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
        G_O = min(gibbs_energies['G_O'])
        
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
            gibbs_energies['dGmax'] = gibbs_energies[['dG1', 'dG2', 'dG3', 'dG4']].idxmax(axis=1)
            gibbs_energies['dGmin'] = gibbs_energies[['dG1', 'dG2', 'dG3', 'dG4']].idxmin(axis=1)
        else:
            gibbs_energies[['OER', 'ORR', 'dGmax', 'dGmin']] = None
            
        gibbs_energies = gibbs_energies.set_index(energies['clean'].index)

        OER = max(dG1, dG2, dG3, dG4) - 1.23
        ORR = 1.23 - min(dG1, dG2, dG3, dG4)
        
        for index in energies['clean'].index:
            spin_cross_over.loc[index, 'clean'] = energies['clean']['spin'].loc[index]
            spin_cross_over.loc[index, 'OH'] = energies['OH']['spin'].loc[index]
            spin_cross_over.loc[index, 'O'] = energies['O']['spin'].loc[index]
        
        gibbs_energies.to_csv(f'{row}_{group}{metal}_gibbs.tsv', sep='\t', float_format='%.2f')
        spin_cross_over.to_csv(f'{row}_{group}{metal}_spin.tsv', sep='\t')
        print(f"Data saved to {row}_{group}{metal}_gibbs.tsv and {row}_{group}{metal}_spin.tsv")
        
        # plotting(gibbs_energies=gibbs_energies, spin_cross_over=spin_cross_over, row=row, group=group, metal=metal,
        #          rxn='OER', rds='dGmax', overpotential=OER, ymin=0.2, ymax=1.4)
        # plotting(gibbs_energies=gibbs_energies, spin_cross_over=spin_cross_over, row=row, group=group, metal=metal,
        #          rxn='ORR', rds='dGmin', overpotential=ORR, ymin=0.2, ymax=1.4)
        # print(f"Figures saved as {row}_{group}{metal}_OER.png and {row}_{group}{metal}_ORR.png")
    
    scaling_relationship['G_'] = scaling_relationship['G_']
    scaling_relationship['G_OH'] = scaling_relationship['G_OH'] + OH_corr
    scaling_relationship['G_O'] = scaling_relationship['G_O'] + O_corr

    scaling_relationship['dG_OH'] = scaling_relationship['G_OH'] - scaling_relationship['G_'] - hydroxide_G
    scaling_relationship['dG_O'] = scaling_relationship['G_O'] - scaling_relationship['G_'] - oxygen_G
    scaling_relationship['dG_OOH'] = scaling_relationship['dG_OH'] + 3.2

    scaling_relationship['dG1'] = scaling_relationship['dG_OH']
    scaling_relationship['dG2'] = scaling_relationship['dG_O'] - scaling_relationship['dG_OH']
    scaling_relationship['dG3'] = scaling_relationship['dG_OOH'] - scaling_relationship['dG_O']
    scaling_relationship['dG4'] = 4.92 - scaling_relationship['dG_OOH']
    
    if scaling_relationship[['dG1', 'dG2', 'dG3', 'dG4']].notna().all().all():
        scaling_relationship['OER'] = scaling_relationship[['dG1', 'dG2', 'dG3', 'dG4']].max(axis=1) - 1.23
        scaling_relationship['ORR'] = 1.23 - scaling_relationship[['dG1', 'dG2', 'dG3', 'dG4']].min(axis=1)
        scaling_relationship['dGmax'] = scaling_relationship[['dG1', 'dG2', 'dG3', 'dG4']].idxmax(axis=1)
        scaling_relationship['dGmin'] = scaling_relationship[['dG1', 'dG2', 'dG3', 'dG4']].idxmin(axis=1)
    else:
        scaling_relationship[['OER', 'ORR', 'dGmax', 'dGmin']] = None
        
    scaling_relationship.to_csv('scaling_relationship.tsv', sep='\t', float_format='%.2f')
    volcano(scaling_relationship, rxn='OER', rds='dGmax', descriptor='dG2', xlabel='O-OH (dG2)', 
            xmin=-2.0, xmax=3.0, ymin=-4.0, ymax=1.0)
    volcano(scaling_relationship, rxn='ORR', rds='dGmin', descriptor='dG1', xlabel='OH (dG1)', 
            xmin=-3.0, xmax=2.0, ymin=-4.0, ymax=1.0)
    
    scaling(scaling_relationship, metals)

def volcano(scaling_relationship, rxn, rds, descriptor, xlabel, xmin, xmax, ymin, ymax):
    x = scaling_relationship[descriptor]
    y = -scaling_relationship[rxn]
    
    if rxn == 'OER':
        y_vals = [-(scaling_relationship[f'dG{i+1}'] - 1.23) for i in range(4)]
    elif rxn == 'ORR':
        y_vals = [-(1.23 - scaling_relationship[f'dG{i+1}']) for i in range(4)]
        
    xx = np.linspace(xmin, xmax, 100)
    plt.figure(figsize=(4, 3), dpi=300)
    
    for i in range(4):
        coeffs = np.polyfit(x, y_vals[i], 1)  # Fit a linear trend
        trendline = np.poly1d(coeffs)
        plt.plot(xx, trendline(xx), label=f'dG{i+1} (trend)', linestyle='-', color=colors[i])
        
    plt.scatter(x, -scaling_relationship[rxn], color='black', s=20, zorder=3)  # Activity vs. Descriptor
    
    # Annotate the points with metal labels
    for xi, yi, metal in zip(x, -scaling_relationship[rxn], metals):
        plt.annotate(f'{metal}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='black')

    # Format axes and limits
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    # Set labels and legend
    plt.xlabel(xlabel)
    plt.ylabel(f'{rxn} activity (-η, eV)')
    plt.legend(labelspacing=0.3)
    
    # Tight layout and save the plot
    plt.tight_layout()
    plt.savefig(f'volcano_{rxn}.png')
    print(f"Figure saved as volcano_{rxn}.png")
    plt.close()

def scaling(scaling_relationship, metals):
    xx = np.linspace(min(scaling_relationship['dG_OH']), max(scaling_relationship['dG_OH']), 100)
    x = scaling_relationship['dG_OH']
    y = scaling_relationship['dG_O']
    plt.figure(figsize=(4.75, fig_height), dpi=300)
    plt.scatter(x, y, c=colors, s=20)
    for xi, yi, metal in zip(x, y, metals):
        plt.annotate(f'{metal}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='black')
    coeffs = np.polyfit(x, y, 1)
    line = np.poly1d(coeffs)
    plt.plot(xx, line(xx), label=r'$\Delta$G$_{\sf O}$ (trend)', linestyle='-', color='black')
    equation = f'y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}'
    plt.text(0.1, 0.8 if coeffs[0] > 0 else 0.1, equation, transform=plt.gca().transAxes, fontsize=10, color='black')
    plt.xlabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize='large')
    plt.ylabel(r'$\Delta$G$_{\sf O}$ (eV)', fontsize='large')
    plt.xlim(-2.5, 1.5)
    plt.ylim(-4.5, 4.5)
    plt.tight_layout()
    plt.savefig('scaling_relationship.png')
    print("Figure saved as scaling_relationship.png")
    plt.close()
        
def plotting(gibbs_energies, spin_cross_over, row, group, metal, rxn, rds, overpotential, ymin, ymax):
    if gibbs_energies.isna().all().all():
        print("Dataframe contains only NaN values, skipping plot.")
        return
    png_filename = f'{row}_{group}{metal}_{rxn}.png'
    marker_size = 0.03
    fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
    plt.xlabel('dz (Å)')
    plt.ylabel(f'{rxn} overpotential (eV)')
    plt.ylim(ymin, ymax)
    plt.yticks(np.arange(ymin, ymax+0.2, 0.2))
    filtered_gibbs_energies = gibbs_energies[rxn].dropna()
    if not filtered_gibbs_energies.empty:
        x = filtered_gibbs_energies.index
        y = filtered_gibbs_energies.values
        try:
            x_new = np.linspace(min(x), max(x), 300)
            spl = make_interp_spline(x, y, k=3) if len(x) > 3 else make_interp_spline(x, y, k=2)
            y_smooth = spl(x_new)
            ax.plot(x_new, y_smooth, color='black', zorder=1)
            ax.scatter(x, y, s=100, color='none', zorder=2)
            for xi, yi in zip(x, y):
                color_ = ms_colors[spin_cross_over.loc[xi, 'clean']]
                color_OH = ms_colors[spin_cross_over.loc[xi, 'OH']]
                color_O = ms_colors[spin_cross_over.loc[xi, 'O']]
                color_OOH = 'white'
                dGrds = gibbs_energies.loc[xi, rds]

                if dGrds == 'dG1':
                    plot_two_color_marker(ax, xi, yi, size=marker_size, color1=color_, color2=color_OH)
                elif dGrds == 'dG2':
                    plot_two_color_marker(ax, xi, yi, size=marker_size, color1=color_OH, color2=color_O)
                elif dGrds == 'dG3':
                    plot_two_color_marker(ax, xi, yi, size=marker_size, color1=color_O, color2=color_OOH)
                elif dGrds == 'dG4':
                    plot_two_color_marker(ax, xi, yi, size=marker_size, color1=color_OOH, color2=color_)
                
                ax.annotate(dGrds, (xi, yi), textcoords="offset points", xytext=(0, 6), ha='center', color='black')
        except ValueError as e:
            print(f"Error while creating spline: {e}")    
    plt.axhline(y=overpotential, color='black', linestyle='--', linewidth=1.0, zorder=0)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tight_layout()
    plt.savefig(png_filename)
    plt.close()
    
def plot_two_color_marker(ax, x, y, size, color1, color2):
    lw = 1.0
    edgecolor = 'black'
    left_wedge = Wedge((x, y), size, 90, 270, facecolor=color1, edgecolor=edgecolor, lw=lw)
    right_wedge = Wedge((x, y), size, 270, 90, facecolor=color2, edgecolor=edgecolor, lw=lw)
    ax.add_patch(left_wedge)
    ax.add_patch(right_wedge)

if __name__ == '__main__':
    main()
