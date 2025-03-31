#!/usr/bin/env python3

import os
import sys
import warnings
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

warnings.filterwarnings('ignore')
png_name = 'nitrogenase'

T = 273.15 + 25
kJmol = 96.485
Jmol = 96.485 * 1E+3
calmol = 23.061
water = 2.4583 # the standard Gibbs free energy of formation of water
kbt = 0.0256 
const = kbt * np.log(10)

# gas
h2 = -6.98952052 # RPBE
zpeh2 = 0.291 # RPBE
cvh2 = 0.091 # RPBE
tsh2 = 0.434 # RPBE

dgh2 = zpeh2 - tsh2 + cvh2
gh2 = h2 + zpeh2 - tsh2 + cvh2
gh = gh2 / 2

h2s = -11.25130575 # RPBE
zpeh2s = 0.400 # RPBE
cvh2s = 0.105 # RPBE
tsh2s = -0.636 # RPBE

dgh2s = zpeh2s - tsh2s + cvh2s
gh2s = h2s + zpeh2s - tsh2s + cvh2s

# solid
gs = -126.00103840/32 # RPBE

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

dgoh = zpeoh + cvoh - tsoh
dgo = zpeo + cvo - tso
dgh = dgoh - dgo
dgs = dgh2s - dgh2 ## caution

print(f"dgh: {dgh:.2f}, dgo: {dgo:.2f}, dgs: {dgs:.2f}")

df = pd.DataFrame()
df.loc['EO', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-166.3055962, 2, 8, 7, 1, 1, 2, 10]
df.loc['E1', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-170.0337744, 2, 9, 7, 1, 1, 2, 10]
df.loc['E2', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-173.8708366, 2, 10, 7, 1, 1, 2, 10]
df.loc['E3', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-177.5720029, 2, 11, 7, 1, 1, 2, 10]
df.loc['E4', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-180.622922, 2, 12, 7, 1, 1, 2, 10]
df.loc['Vac', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-168.4517031, 2, 10, 7, 1, 1, 2, 9] # C2H8Fe7MoNO2S9 + H2

df['comp'] = (
    'X' +
    # 'C' + (df['#C']-2).astype(str) +
    'H' + (df['#H']-8).astype(str) +
    # 'Fe' + (df['#Fe']-7).astype(str) +
    # 'Mo' + (df['#Mo']-1).astype(str) +
    # 'N' + (df['#N']-1).astype(str) +
    # 'O' + (df['#O']-2).astype(str) +
    'S' + (df['#S']-9).astype(str)
)
df['G'] = df['E'] + (0
    # + dgc * (df['#C']-2)
    + dgh * (df['#H']-8)
    # + dgfe * (df['#Fe']-7)
    # + dgmo * (df['#Mo']-1)
    # + dgn * (df['#N']-1)
    # + dgo * (df['#O']-2)
    + dgo * (df['#S']-9) ## caution
) - (0
    # + gc * (df['#C']-2)
    + gh * (df['#H']-8)
    # + gfe * (df['#Fe']-7)
    # + gmo * (df['#Mo']-1)
    # + gn * (df['#N']-1)
    # + go * (df['#O']-2)
    + gs * (df['#S']-9)
)
df['energy'] = df['G'] - df.loc['Vac', 'G'] # - water * df['#O']

df_print = pd.DataFrame()
df_print = df
df_print['G'] = df_print['G'].astype(float).round(2)
df_print['energy'] = df_print['energy'].astype(float).round(2)
print(df_print)

def get_ref_entries():
    ref_entries = []
    
    refs={
        # 'X': 'C2H8Fe7MoNO2S9',
        # 'C': 'C(s)',
        # 'Fe': 'Fe(s)',
        # 'Mo': 'Mo(s)',
        # 'N': 'N2(g)',
        'S': 'S(s)',
        }
    
    for comp, name in refs.items():
        entry = PourbaixEntry(PDEntry(comp, 0.0, name=name))
        ref_entries.append(entry)

    return ref_entries
    
def get_cluster_entries():
    cluster_entries = []
    
    for index, row in df.iterrows():    
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy']))
        cluster_entries.append(entry)
    
    return cluster_entries  

def get_solid_entries():
    solid_entries = []
    solids={
        # 'FeO' : -58.880/calmol,
        # 'Fe3O4': -242.400/calmol,
        # 'Fe2O3': -177.100/calmol,
        # 'Fe2O3': -161.930/calmol,
        # 'Fe(OH)2': -115.570/calmol,
        # 'Fe(OH)3': -166.000/calmol,
        # 'MoO2' : -120.000/calmol,
        # # 'MoO3': -227.000/calmol,
        # 'MoO3': -161.950/calmol,
        # 'FeS': -101.67/kJmol,
        # 'FeS2': -167.36/kJmol,
        # 'FeS2': -171.54/kJmol, # marcasite
        # 'Mo2S3': -407.10/kJmol,
        # 'MoS2': -276.14/kJmol,
        }
    
    for solid, energy in solids.items():
        entry = PourbaixEntry(PDEntry(solid, energy))
        solid_entries.append(entry)

    return solid_entries
          
def get_gas_entries():
    gas_entries = []
    gases={
        # # 'CO': -110.53/kJmol,
        # # 'CO2': -393.52/kJmol,
        # # 'CH4': -74.87/kJmol,
        # # 'SO2': -296.84/kJmol-T*248.223/Jmol,
        # # 'SH2': -20.6/kJmol-T*205.81/Jmol,
        # # 'NO': 90.29/kJmol,
        # # 'NO2': 33.10/kJmol,
        # # 'NH3': -45.94/kJmol, # All from NIST Chem WebBook
        'H2S': -7.892/calmol,
        'SO': 12.780/calmol,
        'SO2': -71.790/calmol,
        'SO3': -88.520/calmol,
        }
    
    for gas, energy in gases.items():
        entry = PourbaixEntry(PDEntry(gas, energy))
        gas_entries.append(entry)

    return gas_entries
          
def get_ion_entries(concentration=1e-6):
    ion_entries = []
    ions={
        # 'Fe++': -20.300/calmol,
        # 'HFeO2-': -90.627/calmol,
        # 'Fe+++': -2.530/calmol,
        # 'FeOH++': -55.910/calmol,
        # 'Fe(OH)2+': -106.200/calmol,
        # # 'FeO4-': -111685/calmol,
        # 'Mo+++': -13.800/calmol,
        # # 'HMoO4-': -213.600/calmol,
        # 'MoO4--': -205.420/calmol,
        'H2S': -6.540/calmol,
        'HS-': 3.010/calmol,
        'S--': 21.958/calmol,
        'S2--': 19.749/calmol,
        'S3--': 17.968/calmol,
        'S4--': 16.615/calmol,
        'S5--': 15.689/calmol,
        'H2S2O3': -129.900/calmol,
        'HS2O3-': -129.500/calmol,
        'S2O3--': -127.200/calmol,
        'S5O6--': -228.500/calmol,
        'S4O6--': -244.300/calmol,
        'HS2O4-': -141.408/calmol,
        'S2O4--': -138.000/calmol,
        'S3O6--': -229.000/calmol,
        'H2SO3': -128.690/calmol,
        'HSO3-': -126.000/calmol,
        'SO3--': -116.100/calmol,
        'S2O6--': -231.000/calmol,
        'H2SO4': -177.340/calmol,
        'HSO4-': -179.940/calmol,
        'SO4--': -177.340/calmol,
        'S2O8--': -262.000/calmol,
        }

    # updated_ions = {ion: energy - 0.05917 * log10(concentration) for ion, energy in ions.items()}
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy))
        ion_entries.append(entry)

    return ion_entries
    
def plot_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)
    stable_entries = pourbaix.stable_entries
    
    fig, ax = plt.subplots(figsize=(6, 5))    
    plotter.get_pourbaix_plot(limits=[[0, 14], [-2, 1]], label_domains=False,
                              show_water_lines=False, show_neutral_axes=False, ax=ax)
    
    for line in ax.lines:
        line.set_linewidth(1.0)
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_color('black')
        text.set_fontweight('bold')
        
    pH2 = np.arange(0, 14.01, 0.01)
    plt.plot(pH2, 0.05 - pH2 * const, color='blue', lw=1.5, dashes=(3, 1))
    plt.plot(pH2, -0.7 - pH2 * const, color='black', lw=1.5, dashes=(3, 1))
    # ax.text(7, -0.5, r'E°$_{N_{2}/NH_{3}}$ = 0.05 V vs RHE', 
    #         color='blue', rotation=-14, fontsize=14, ha='center')
    ax.text(2.5, -0.45, r'E°$_{N_{2}/NH_{3}}$ = 0.05 V vs RHE', 
            color='blue', rotation=-14, fontsize=14)
    
    name_mapping1 = {
        'H2SO4(aq) + XH2(s)': 'XH2(s) + H2SO4(aq)',
        'HSO4[-1] + XH2(s)': 'XH2(s) + HSO4[-1]',
        'SO4[-1] + XH2(s)': 'XH2(s) + SO4[-1]',
        'SO4[-2] + XH2(s)': 'XH2(s) + SO4[-2]',
        'H2S(aq) + XH2(s)': 'XH2(s) + H2S(aq)',
        'HS[-1] + XH2(s)': 'XH2(s) + HS[-1]',
        'S[-2] + XH2(s)': 'XH2(s) + S[-2]',
    }
    
    name_mapping2 = {
        # 'X(s)': r'E$_{vv}$',
        'XS(s)': r'E$_{0}$', 
        'XHS(s)': r'E$_{1}$',
        'XH2S(s)': r'E$_{2}$',
        'XH3S(s)': r'E$_{3}$',
        'XH4S(s)': r'E$_{4}$',
        'XH2(s)': r'E$_{v}$',
        'S[-2]': 'S²⁻',
        'HS[-1]': 'HS⁻',
        'H2S(aq)': 'H₂S(aq)',
        ' SO4[-1]': ' S₂O₈²⁻',
        'SO4[-2]': 'SO₄²⁻',
        'HSO4[-1]': 'HSO₄⁻',
        'H2SO4(aq)': 'H₂SO₄(aq)',
    }

    for text in ax.texts:
        old_name = text.get_text()
        new_name = old_name
        for old_part, new_part in name_mapping1.items():
            if old_part in new_name:
                new_name = new_name.replace(old_part, new_part)
        text.set_text(new_name)
        
    for text in ax.texts:
        old_name = text.get_text()
        new_name = old_name
        for old_part, new_part in name_mapping2.items():
            if old_part in new_name:
                new_name = new_name.replace(old_part, new_part)
        text.set_text(new_name)
        
    sulfur_names = [
        'XS(s)',
        'XHS(s)',
        'XH2S(s)',
        'XH3S(s)',
        'XH4S(s)',
    ]
    
    vacancy_names = [
        'HSO4[-1]',
        'SO4[-1]',
        'HS[-1]',
        'H2SO4(aq)',
        'SO4[-2]',
        'S[-2]',
        'H2S(aq)',
    ]
        
    sulfur_entries = [entry for entry in stable_entries if 'XH2(s)' not in entry.name]
    vacancy_entries = [entry for entry in stable_entries if 'XH2(s)' in entry.name]
    sulfur_colors = [plt.cm.Greys(i) for i in np.linspace(0.05, 0.35, len(sulfur_names))]
    vacancy_colors = [plt.cm.Blues(i) for i in np.linspace(0.05, 0.35, len(vacancy_names))] # YlOrBr

    for i, entry in enumerate(sulfur_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        color = sulfur_colors[-sulfur_names.index(entry.name)]
        ax.fill(x, y, color=color)
        
    for i, entry in enumerate(vacancy_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        ion_part = entry.name.replace(' + ', '').replace('XH2(s)', '')
        color = vacancy_colors[vacancy_names.index(ion_part)]
        ax.fill(x, y, color=color)
        
    if 'cluster' in png_name:
        ax.text(13, -0.4, r"$\mathbf{E_{0}}$", fontsize=14, ha='center')
        ax.text(13, -0.9, r"$\mathbf{E_{2}}$", fontsize=14, ha='center')
        ax.text(13, -1.25, r"$\mathbf{E_{3}}$", fontsize=14, ha='center')
        ax.text(13, -1.85, r"$\mathbf{E_{4}}$", fontsize=14, ha='center')
    elif 'bulk' in png_name:
        ax.text(12, -0.50, r"$\mathbf{E_{0}}$", fontsize=14, ha='center')
        ax.text(12, -0.80, r"$\mathbf{E_{2}}$", fontsize=14, ha='center')
        ax.text(12, -1.00, r"$\mathbf{E_{3}}$", fontsize=14, ha='center')
        ax.text(0.7, 0.8, r"$\mathbf{E_{v} + HSO_{4}^{-}}$", fontsize=14)
        ax.text(9.0, 0.5, r"$\mathbf{E_{v} + SO_{4}^{2-}}$", fontsize=14)
        ax.text(3.5, -1.6, r"$\mathbf{E_{v} + H_{2}S(aq)}$", fontsize=14, ha='center')
        ax.text(9.0, -1.5, r"$\mathbf{E_{v} + HS^{-}}$", fontsize=14)
        ax.text(13.8, -1.9, r"$\mathbf{E_{v} + S^{2-}}$", fontsize=14, ha='right')
        
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
        
    plt.tight_layout()
    plt.savefig(png_name, dpi=300, bbox_inches='tight')
    plt.show()

def main():
    print('\n################## Reference Entries ##########################################\n')
    ref_entries = get_ref_entries()
    for entry in ref_entries:
        print(entry)
    
    print('\n################## cluster Entries ##########################################\n')
    cluster_entries = get_cluster_entries()
    for entry in cluster_entries:
        print(entry)

    print('\n################## Solid Entries ##########################################\n')
    solid_entries = get_solid_entries()
    for entry in solid_entries:
        print(entry)

    print('\n################## Gas Entries ##########################################\n')
    gas_entries = get_gas_entries()
    for entry in gas_entries:
        print(entry)
        
    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = ref_entries + cluster_entries + solid_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))
    
    all_entries = ref_entries + cluster_entries
    plot_pourbaix(all_entries, f'{png_name}_cluster4.png')

    # for i in range(7):
    #     ion_entries = get_ion_entries(concentration=10**(-i))
    #     all_entries = solid_entries + gas_entries + ion_entries
    #     # for entry in all_entries:
    #     #     print(entry.conc_term)
    #     plot_pourbaix(all_entries, f'{png_name}_ion{i}.png')
    
    # plot_pourbaix(solid_entries, f'{png_name}_solid.png')
    # all_entries = solid_entries + ion_entries
    # plot_pourbaix(all_entries, f'{png_name}_exp.png')

    # plot_pourbaix(mpr1_entries, f'{png_name}_mpr1.png')
    # plot_pourbaix(mpr2_entries, f'{png_name}_mpr2.png')
    
    all_entries = ref_entries + cluster_entries + solid_entries + gas_entries + ion_entries # + mpr_entries
    plot_pourbaix(all_entries, f'{png_name}_bulk4.png')

if __name__ == "__main__":
    main()
