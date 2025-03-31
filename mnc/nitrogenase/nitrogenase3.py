#!/usr/bin/env python3

import os
import sys
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

warnings.filterwarnings('ignore')
png_name = 'iron_sulfur'

# API_KEY = os.getenv('MAPI_KEY')
# if not API_KEY:
#     sys.exit("Error: MAPI_KEY environment variable not set.")
# mpr = MPRester(API_KEY)
# mpr_entries = mpr.get_pourbaix_entries(['Fe', 'S', 'Mo', 'C'])

# mpr1_entries = []
# mpr2_entries = []

# for entry in mpr_entries:
#     if 'ion' in entry.entry_id and entry.npH - entry.nPhi > 0:
#         mpr1_entries.append(entry)
#     else:
#         mpr2_entries.append(entry) 

T = 273.15 + 25
kJmol = 96.485
Jmol = 96.485 * 1E+3
calmol = 23.061
water = 2.4583 # the standard Gibbs free energy of formation of water

# gas
h2 = -6.98952052
h2o = -14.1655297

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

zpeh2o = 0.560
cvh2o = 0.103
tsh2o = 0.675

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2

# n2 = -16.64503942
# zpen2 = 0.098
# tsn2 = 0.592
# gn2 = n2 + zpen2 - tsn2

gh = gh2 / 2
go = gh2o - gh2
# gn = gn2 / 2

# solid
# gc = -.37429454E+02/4
# gfe = -5.041720865
# gmo = -10.94947049
gs = -126.00103840/32

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

df = pd.DataFrame()
df.loc['EO', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-166.3055962, 2, 8, 7, 1, 1, 2, 10]
df.loc['E1', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-170.0337744, 2, 9, 7, 1, 1, 2, 10]
df.loc['E2', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-173.8708366, 2, 10, 7, 1, 1, 2, 10]
df.loc['E3', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-177.5720029, 2, 11, 7, 1, 1, 2, 10]
df.loc['E4', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-180.622922, 2, 12, 7, 1, 1, 2, 10]
df.loc['Vac', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-168.4517031, 2, 10, 7, 1, 1, 2, 9] # C2H8Fe7MoNO2S9 + H2
reference = -168.4517031 + 2 * dgh - gh2 # C2H8Fe7MoNO2S9

df['comp'] = (
    'X' +
    # 'C' + df['#C'].astype(str) +
    'H' + (df['#H']-8).astype(str) +
    # 'Fe' + df['#Fe'].astype(str) +
    # 'Mo' + df['#Mo'].astype(str) +
    # 'N' + df['#N'].astype(str) +
    # 'O' + df['#O'].astype(str) +
    'S' + (df['#S']-9).astype(str)
)
df['G'] = df['E'] + (
    dgh * (df['#H']-8)
    # + dgo * df['#O']
) - (
    # gc * df['#C'] 
    + gh * (df['#H']-8)
    # + gfe * df['#Fe']
    # + gmo * df['#Mo']
    # + gn * df['#N']
    # + go * df['#O']
    + gs * (df['#S']-9)
) 
df['energy'] = df['G'] - reference # - water * df['#O']

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
    # solids={
    #     'FeO' : -58.880/calmol,
    #     'Fe3O4': -242.400/calmol,
    #     'Fe2O3': -177.100/calmol,
    #     'Fe2O3': -161.930/calmol,
    #     'Fe(OH)2': -115.570/calmol,
    #     'Fe(OH)3': -166.000/calmol,
    #     'MoO2' : -120.000/calmol,
    #     # 'MoO3': -227.000/calmol,
    #     'MoO3': -161.950/calmol,
    #     'FeS': -101.67/kJmol,
    #     'FeS2': -167.36/kJmol,
    #     'FeS2': -171.54/kJmol, # marcasite
    #     'Mo2S3': -407.10/kJmol,
    #     'MoS2': -276.14/kJmol,
    #     }
    
    # for solid, energy in solids.items():
    #     entry = PourbaixEntry(PDEntry(solid, energy))
    #     solid_entries.append(entry)

    return solid_entries
          
def get_gas_entries():
    gas_entries = []
    gases={
        # 'CO': -110.53/kJmol,
        # 'CO2': -393.52/kJmol,
        # 'CH4': -74.87/kJmol,
        'SO2': -296.84/kJmol-T*248.223/Jmol,
        'SH2': -20.6/kJmol-T*205.81/Jmol,
        # 'NO': 90.29/kJmol,
        # 'NO2': 33.10/kJmol,
        # 'NH3': -45.94/kJmol, # All from NIST Chem WebBook
        }
    
    for gas, energy in gases.items():
        entry = PourbaixEntry(PDEntry(gas, energy))
        gas_entries.append(entry)

    return gas_entries
          
def get_ion_entries():
    ion_entries = []
    # ions={
    #     'Fe++': -20.300/calmol,
    #     'HFeO2-': -90.627/calmol,
    #     'Fe+++': -2.530/calmol,
    #     'FeOH++': -55.910/calmol,
    #     'Fe(OH)2+': -106.200/calmol,
    #     # 'FeO4-': -111685/calmol,
    #     'Mo+++': -13.800/calmol,
    #     # 'HMoO4-': -213.600/calmol,
    #     'MoO4--': -205.420/calmol,
    #     }
    
    # for ion, energy in ions.items():
    #     comp = Ion.from_formula(ion)
    #     entry = PourbaixEntry(IonEntry(comp, energy), concentration=1e-6)
    #     ion_entries.append(entry)

    return ion_entries
    
def plot_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    # ax = plotter.get_pourbaix_plot(limits=[[0, 14], [-2, 1.5]])
    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-4, 4]])
    
    for line in ax.lines:
        line.set_linewidth(1.0)
    for text in ax.texts:
        text.set_fontsize(14)
    
    name_mapping1 = {
        # 'X(s)': r'E$_{vv}$',, 
        # 'XS(s)': r'E$_{0}$', 
        # 'XHS(s)': r'E$_{1}$',
        # 'XH2S(s)': r'E$_{2}$',
        # 'XH3S(s)': r'E$_{3}$',
        # 'XH4S(s)': r'E$_{4}$',
        # 'XH2(s)': r'E$_{v}$',
        # 'SO2(s)': 'SO₂(s)',
        # 'H2S(s)': 'H₂S(s)',
        'SO2(s) + X(s)': 'X(s) + SO2(s)',
        'H2S(s) + XH2(s)': 'XH2(s) + H2S(s)',
    }
    
    name_mapping2 = {
        #₁₂₃₄₅₆₇₈₉₀
        'X(s)': 'Fe₇MoC₂NO₂H₈S₉',  
        'XS(s)': 'Fe₇MoC₂NO₂H₈S₁₀', 
        'XHS(s)': 'Fe₇MoC₂NO₂H₉S₁₀', 
        'XH2S(s)': 'Fe₇MoC₂NO₂H₁₀S₁₀', 
        'XH3S(s)': 'Fe₇MoC₂NO₂H₁₁S₁₀', 
        'XH4S(s)': 'Fe₇MoC₂NO₂H₁₂S₁₀', 
        'XH2(s)': 'Fe₇MoC₂NO₂H₁₀S₉',
        '2': '₂',
        '3': '₃',
        '4': '₄',
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
        
    # if 'cluster' in png_name:
    #     omit_parts = ['X(s)', 'Fe(s)', 'Mo(s)', 'N2(s)', 'C(s)' , 'S(s)', ' ', '+']
    #     for text in ax.texts:
    #         old_name = text.get_text()
    #         new_name = old_name
    #         for old_part in omit_parts:
    #             if old_part in old_name:
    #                 new_name = new_name.replace(old_part, '')
    #             text.set_text(new_name)
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    fig = ax.figure
    fig.set_size_inches((10, 8))
    
    plt.savefig(png_name, dpi=300, bbox_inches='tight') #, transparent=True)
    # plt.close()
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
    plot_pourbaix(all_entries, f'{png_name}_cluster3.png')
    
    # plot_pourbaix(solid_entries, f'{png_name}_solid.png')
    # plot_pourbaix(ion_entries, f'{png_name}_ion.png')
    # all_entries = solid_entries + ion_entries
    # plot_pourbaix(all_entries, f'{png_name}_exp.png')

    # plot_pourbaix(mpr1_entries, f'{png_name}_mpr1.png')
    # plot_pourbaix(mpr2_entries, f'{png_name}_mpr2.png')
    
    all_entries = ref_entries + cluster_entries + solid_entries + gas_entries + ion_entries
    plot_pourbaix(all_entries, f'{png_name}_bulk3.png')

if __name__ == "__main__":
    main()
