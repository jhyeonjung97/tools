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
import argparse

warnings.filterwarnings('ignore')
png_name = 'protonation'
tag = '_all'

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
tsh2s_aq = 0.93 # RPBE, aq
tsh2s_gas = 0.636 # RPBE, aq

dgh2s_aq = zpeh2s - tsh2s_aq + cvh2s
dgh2s_gas = zpeh2s - tsh2s_gas + cvh2s
gh2s_aq = h2s + zpeh2s - tsh2s_aq + cvh2s
gh2s_gas = h2s + zpeh2s - tsh2s_gas + cvh2s

# solid
s = -126.00103840/32 # RPBE
gs = s
fh2s_aq = gh2s_aq - gh2 - gs
fh2s_gas = gh2s_gas - gh2 - gs

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

dgoh = zpeoh + cvoh - tsoh
dgo = zpeo + cvo - tso
dgh = 0.217

df = pd.DataFrame()
df.loc['E₀',     ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-412.4041388, 0, 0, 0, 0, 0, 0,  0]
df.loc['E₀*H',   ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-416.8018053, 0, 1, 0, 0, 0, 0,  0]
df.loc['E₀*H₂',  ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-419.8147836, 0, 2, 0, 0, 0, 0,  0]
df.loc['E₁',     ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-416.1486068, 0, 1, 0, 0, 0, 0,  0]
df.loc['E₁*H',   ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-420.5210881, 0, 2, 0, 0, 0, 0,  0]
df.loc['E₁*H₂',  ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-423.4465657, 0, 3, 0, 0, 0, 0,  0]
df.loc['E₂',     ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-420.0688973, 0, 2, 0, 0, 0, 0,  0]
df.loc['E₂*H',   ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-424.3006724, 0, 3, 0, 0, 0, 0,  0]
df.loc['E₂*H₂',  ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-427.0052327, 0, 4, 0, 0, 0, 0,  0]
df.loc['Vac',    ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-414.4369611, 0, 2, 0, 0, 0, 0, -1]
df.loc['Vac*H',  ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-418.5229228, 0, 3, 0, 0, 0, 0, -1]
df.loc['Vac*H₂', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-421.2139338, 0, 4, 0, 0, 0, 0, -1]
df.loc['E₃',     ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-423.8152177, 0, 3, 0, 0, 0, 0,  0]
df.loc['E₃*H',   ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-427.9188670, 0, 4, 0, 0, 0, 0,  0]
df.loc['E₃*H₂',  ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-430.6141604, 0, 5, 0, 0, 0, 0,  0]

df['comp'] = (
    'X'
    + 'H' + (df['#H']).astype(str)
    + 'S' + (df['#S']+1).astype(str)
)
df['G'] = df['E'] + dgh * df['#H'] - gh * df['#H'] - gs * df['#S']
df['energy'] = df['G'] - df.loc['E₀', 'G'] # - water * df['#O']
print(df)

def get_ref_entries():
    ref_entries = []
    entry = PourbaixEntry(PDEntry('S', 0.0))
    entry = PourbaixEntry(PDEntry('H2S', fh2s_aq))
    ref_entries.append(entry)
    
    return ref_entries
    
def get_cluster_entries():
    cluster_entries = []
    
    for index, row in df.iterrows():    
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy']))
        cluster_entries.append(entry)
    
    return cluster_entries  
          
def get_gas_entries():
    gas_entries = []
    gases={
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
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy))
        ion_entries.append(entry)

    return ion_entries
    
def plot_2d_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)
    stable_entries = pourbaix.stable_entries
    
    fig, ax = plt.subplots(figsize=(6, 5))    
    plotter.get_pourbaix_plot(limits=[[0, 14], [-2, 1]], label_domains=True,
                              show_water_lines=False, show_neutral_axes=False, ax=ax)
    
    for line in ax.lines:
        line.set_linewidth(1.0)
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_color('black')
        
    pH2 = np.arange(0, 14.01, 0.01)
    plt.plot(pH2, 0.05 - pH2 * const, color='blue', lw=1.5, dashes=(3, 1))
    plt.plot(pH2, -0.7 - pH2 * const, color='black', lw=1.5, dashes=(3, 1))
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
        'H2S(s) + XH3(s)': 'XH3(s) + H2S(s)',
        'H2S(s) + XH4(s)': 'XH4(s) + H2S(s)',
        'HS[-1] + XH3(s)': 'XH3(s) + HS[-1]',
        'HS[-1] + XH4(s)': 'XH4(s) + HS[-1]',
        'S[-2] + XH3(s)': 'XH3(s) + S[-2]',
        'S[-2] + XH4(s)': 'XH4(s) + S[-2]',
    }

    name_mapping2 = {
        ' S[-2]': ' S²⁻',
        ' HS[-1]': ' HS⁻',
        ' H2S(s)': ' H₂S(aq)',
        ' SO4[-1]': ' S₂O₈²⁻',
        ' SO4[-2]': ' SO₄²⁻',
        ' HSO4[-1]': ' HSO₄⁻',
        ' H2SO4(aq)': ' H₂SO₄(aq)',
    }

    count_mapping = {}
    for _, row in df.iterrows():
        comp = row['comp']
        count_mapping[comp] = count_mapping.get(comp, 0) + 1
    
    for comp, count in count_mapping.items():
        if count == 1:
            matching_index = df[df['comp'] == comp].index[0]
            name_mapping2[comp] = matching_index
        else:
            matching_rows = df[df['comp'] == comp]
            min_energy_index = matching_rows.loc[matching_rows['E'].idxmin()].name
            name_mapping2[comp] = min_energy_index

    name_mapping2 = {k.replace('S0', ''): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('S1', 'S'): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('H0', ''): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('H1', 'H'): v for k, v in name_mapping2.items()}
    
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
        'XH5S(s)',
        'XH2(s)',
        'XH3(s)',
        'XH4(s)',
    ]
    
    vacancy_names2 = [
        'SO4[-1]',
        'SO4[-2]',
        'HSO4[-1]',
        'H2SO4(aq)',
    ]

    vacancy_names3 = [
        'H2S(s)',
        'HS[-1]',
        'S[-2]',
    ]

    vacancy_names4 = [
        'H2S(s)',
        'HS[-1]',
        'S[-2]',
    ]
        
    sulfur_entries = [entry for entry in stable_entries if 'XH2(s)' not in entry.name and 'XH3(s)' not in entry.name and 'XH4(s)' not in entry.name]
    vacancy_entries = [entry for entry in stable_entries if 'XH2(s)' in entry.name or 'XH3(s)' in entry.name or 'XH4(s)' in entry.name]
    sulfur_colors = [plt.cm.Greys(i) for i in np.linspace(0.05, 0.35, len(sulfur_names))]
    vacancy_colors2 = [plt.cm.Greens(i) for i in np.linspace(0.15, 0.45, len(vacancy_names2))] # YlOrBr
    vacancy_colors3 = [plt.cm.Blues(i) for i in np.linspace(0.35, 0.45, len(vacancy_names3))] # YlOrBr
    vacancy_colors4 = [plt.cm.Oranges(i) for i in np.linspace(0.35, 0.45, len(vacancy_names4))] # YlOrBr   

    for i, entry in enumerate(sulfur_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        color = sulfur_colors[-sulfur_names.index(entry.name)]
        ax.fill(x, y, color=color)
    
    for i, entry in enumerate(vacancy_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        ion_part = entry.name.replace(' + ', '').replace('XH2(s)', '').replace('XH3(s)', '').replace('XH4(s)', '')
        if 'XH2(s)' in entry.name:
            color = vacancy_colors2[vacancy_names2.index(ion_part)]
        elif 'XH3(s)' in entry.name:
            color = vacancy_colors3[vacancy_names3.index(ion_part)]
        elif 'XH4(s)' in entry.name:
            color = vacancy_colors4[vacancy_names4.index(ion_part)]
        ax.fill(x, y, color=color)
        
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
        
    plt.tight_layout()
    plt.savefig(png_name, dpi=300, bbox_inches='tight')
    print(f'{png_name} saved')
    plt.show()

def plot_1d_pourbaix(entries, png_name, ph=0):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    stable_entries = pourbaix.stable_entries
    
    stable_entries = sorted(stable_entries, key=lambda x: x.name)
    
    print('\n################## Final Entries ##########################################\n')
    for entry in stable_entries:
        print(entry, entry.composition)
    
    fig = plt.figure(figsize=(7, 5), dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([-2.0, 1.0, -6, 2])
    ax.set_xlabel(r'SHE (V)')
    ax.set_ylabel(r'$\Delta$G (eV)')
    
    xx = np.arange(-2.00, 1.05, 0.01)
    
    n_entries = len(stable_entries)
    cmap = plt.cm.tab20
    colors = [cmap(i % 20) for i in range(n_entries)]
    
    name_mapping1 = {
        'H2SO4(aq) + XH2(s)': 'XH2(s) + H2SO4(aq)',
        'HSO4[-1] + XH2(s)': 'XH2(s) + HSO4[-1]',
        'SO4[-1] + XH2(s)': 'XH2(s) + SO4[-1]',
        'SO4[-2] + XH2(s)': 'XH2(s) + SO4[-2]',
        'H2S(aq) + XH2(s)': 'XH2(s) + H2S(aq)',
        'HS[-1] + XH2(s)': 'XH2(s) + HS[-1]',
        'S[-2] + XH2(s)': 'XH2(s) + S[-2]',
        'H2S(s) + XH3(s)': 'XH3(s) + H2S(s)',
        'H2S(s) + XH4(s)': 'XH4(s) + H2S(s)',
        'HS[-1] + XH3(s)': 'XH3(s) + HS[-1]',
        'HS[-1] + XH4(s)': 'XH4(s) + HS[-1]',
        'S[-2] + XH3(s)': 'XH3(s) + S[-2]',
        'S[-2] + XH4(s)': 'XH4(s) + S[-2]',
    }

    name_mapping2 = {
        ' S[-2]': ' S²⁻',
        ' HS[-1]': ' HS⁻',
        ' H2S(s)': ' H₂S(aq)',
        ' SO4[-1]': ' S₂O₈²⁻',
        ' SO4[-2]': ' SO₄²⁻',
        ' HSO4[-1]': ' HSO₄⁻',
        ' H2SO4(aq)': ' H₂SO₄(aq)',
    }

    count_mapping = {}
    for _, row in df.iterrows():
        comp = row['comp']
        count_mapping[comp] = count_mapping.get(comp, 0) + 1
    
    for comp, count in count_mapping.items():
        if count == 1:
            matching_index = df[df['comp'] == comp].index[0]
            name_mapping2[comp] = matching_index
        else:
            matching_rows = df[df['comp'] == comp]
            min_energy_index = matching_rows.loc[matching_rows['E'].idxmin()].name
            name_mapping2[comp] = min_energy_index

    name_mapping2 = {k.replace('S0', ''): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('S1', 'S'): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('H0', ''): v for k, v in name_mapping2.items()}
    name_mapping2 = {k.replace('H1', 'H'): v for k, v in name_mapping2.items()}
    
    for idx, entry in enumerate(stable_entries):
        color = colors[idx]
        name = entry.name        
        dG0 = entry.energy
        comp = entry.composition

        for old_part, new_part in name_mapping1.items():
            if old_part in name:
                name = name.replace(old_part, new_part)
        for old_part, new_part in name_mapping2.items():
            if old_part in name:
                name = name.replace(old_part, new_part)

        dG = dG0 + entry.nPhi * xx + const * entry.npH * ph # SHE
        # dG = dG0 + entry.nPhi * (xx - const * ph) + const * entry.npH * ph # RHE
        print(entry.name, entry.energy, ph, dG[-1])
        ax.plot(xx, dG, '-', lw=1, c=color, label=name)
    
    plt.xlim(-2.0, 1.0)
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02),
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    
    plt.savefig(png_name, bbox_inches='tight')
    print(f'{png_name} saved')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot Pourbaix diagram')
    parser.add_argument('-ph', type=float, default=0.0, help='pH value for 1D plot (default: 7.0)')
    args = parser.parse_args()

    print('\n################## Reference Entries ##########################################\n')
    ref_entries = get_ref_entries()
    for entry in ref_entries:
        print(entry)
    
    print('\n################## cluster Entries ##########################################\n')
    cluster_entries = get_cluster_entries()
    for entry in cluster_entries:
        print(entry)

    print('\n################## Gas Entries ##########################################\n')
    gas_entries = get_gas_entries()
    for entry in gas_entries:
        print(entry)
        
    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = ref_entries + cluster_entries + gas_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))
    plot_1d_pourbaix(all_entries, f'{png_name}_1D{tag}.png', args.ph)
    plot_2d_pourbaix(all_entries, f'{png_name}_2D{tag}.png')
    
if __name__ == "__main__":
    main()
