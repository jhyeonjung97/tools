#!/usr/bin/env python3

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

warnings.filterwarnings('ignore')

filename = '2Co'
tag = ''

kbt = 0.0256 
const = kbt * np.log(10)

kJmol = 96.485
calmol = 23.061
water = 2.4583 # the standard Gibbs free energy of formation of water

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2

n2 = -16.64503942
zpen2 = 0.098
tsn2 = 0.592
gn2 = n2 + zpen2 - tsn2

# solid
gc = -.37429454E+02/4

N4C26 = -.27195317E+03
H2N4C26 = -.28183609E+03

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo

metal_path = '~/bin/tools/mnc/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
gm = metal_df.loc['Co', 'energy']

df = pd.read_csv(f'{filename}_energies.csv', delimiter=',', index_col=0)
df_oer = pd.read_csv(f'{filename}_oer.tsv', delimiter='\t', index_col=0)
df_orr = pd.read_csv(f'{filename}_orr.tsv', delimiter='\t', index_col=0)

df['name'] = 'CoNC(' + df.index.str.upper() + ')'
df['comp'] = 'CoX' + df.index.str.upper().str.replace("-", "")
df['comp'] = df['comp'].str.replace('CoXVAC', 'H2X')
df['name'] = df['name'].str.replace('CoNC(VAC)', 'Co⁺²+H₂NC', regex=False)
df['comp'] = df['comp'].str.replace('CoXCLEAN', 'CoX')
df['name'] = df['name'].str.replace('CoNC(CLEAN)', 'CoNC(clean)')
df['comp'] = df['comp'].str.replace('CoXMH', 'CoXH')
df['comp'] = df['comp'].str.replace('CoXNH', 'CoXH')

# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - 2 * gn2 - 26 * gc - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - N4C26 - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
df['energy'] = df['dG'] + df.loc['clean', 'G'] + gh2 - gm - H2N4C26 - 2 * dgh - water * (df['#O'] + df['#OH'] + df['#OOH']*2)

df = df.drop(index='vac')
df = df.drop(index='oo')
df = df.drop(index='ooh')
df = df.drop(index='o-ooh')
df = df.drop(index='ooh-o')
df = df.drop(index='oh-ooh')
df = df.drop(index='ooh-oh')
df = df.drop(index='ooh-ooh')
# df = df.drop(index='o-o')
# df = df.drop(index='oh-o')
# df = df.drop(index='oh-oh')
df = df.drop(index='ooho')
df = df.drop(index='oohoh')
df = df.drop(index='oohooh')
df = df.dropna()
print(df)

def get_ref_entries():
    ref_entries = []
    
    refs={
        'Co': 'Co(s)',
        # 'N2': 'N2(g)',
        # 'C': 'C(s)',
        # 'X': 'N4C26',
        'H2X': 'H2NC(vac)',
        }
    
    for comp, name in refs.items():
        entry = PourbaixEntry(PDEntry(comp, 0.0, name=name))
        ref_entries.append(entry)
    
    return ref_entries
    
def get_sac_entries():
    sac_entries = []
    
    for index, row in df.iterrows():    
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy'], entry_id=row['name']))
        sac_entries.append(entry)
        
    energy = H2N4C26 + 2 * dgh - gh2 - N4C26
    entry = PourbaixEntry(ComputedEntry('X', -energy, entry_id='NC(vac)'))
    sac_entries.append(entry)
    
    return sac_entries

def get_solid_entries():
    solid_entries = []
    solids={
        'Co': 0,
        'CoO': -52.310/calmol,
        'CoO': -49.000/calmol,
        'Co3O4': -167.835/calmol,
        # 'Co2O3': -115.130/calmol,
        'CoO2': -51.840/calmol,
        'Co(OH)2': -109.999/calmol,
        'Co(OH)3': -142.600/calmol
        }
    
    for solid, energy in solids.items():
        entry = PourbaixEntry(ComputedEntry(solid, energy))
        solid_entries.append(entry)

    return solid_entries

def get_ion_entries():
    ion_entries = []
    ions={
        'Co++': -12.800/calmol,
        'HCoO2-': -82.970/calmol,
        'Co+++': 28.900/calmol,
        }
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy), concentration=1e-6)
        ion_entries.append(entry)

    return ion_entries
    
def get_equation(a, b, c, d, df_rxn, n):
    eq = fr"S$_{{{a}}}$$\rightarrow$S$_{{{b}}}$$\rightarrow$S$_{{{c}}}$$\rightarrow$S$_{{{d}}}$: " + f"{df_rxn['overP'][n]:.2f} eV"
    return eq
    
def plot_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    fig, ax = plt.subplots(figsize=(7, 5))    
    plotter.get_pourbaix_plot(limits=[[0, 14], [-1, 3]], label_domains=False, label_fontsize=14,
                              show_water_lines=False, show_neutral_axes=False, ax=ax)
    stable_entries = pourbaix.stable_entries

    for line in ax.lines:
        line.set_linewidth(0.5)
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_color('black')
        text.set_fontweight('bold')

    if 'bulk' in png_name:
        pH2 = np.arange(0, 14.01, 0.01)
        plt.plot(pH2, 1.23 - pH2 * const, color='black', lw=2.0) #, dashes=(5, 2))
        plt.plot(pH2, df_oer['onsetP'][0] - pH2 * const, color='red', lw=2.0)
        plt.plot(pH2, df_oer['onsetP'][2] - pH2 * const, color='red', lw=2.0, dashes=(5, 2))
        plt.plot(pH2, df_orr['onsetP'][0] - pH2 * const, color='blue', lw=2.0)
        
    vac_entries = [entry for entry in stable_entries if 'XCo' not in entry.name]
    sac_entries = [entry for entry in stable_entries if 'XCo' in entry.name]
    vac_colors = [plt.cm.Greys(i) for i in np.linspace(0.1, 0.3, len(vac_entries))]
    sac_colors = [plt.cm.Blues(i) for i in np.linspace(0.1, 0.3, len(sac_entries))]

    vac_mapping = {
        'XH2(s) + Co(s)': 0,
        'XH2(s) + Co[+2]': 1,
        'X(s) + Co[+2]': 2,
        'X(s) + Co[+3]': 3,
        'Co(s) + XH2(s)': 0,
        'Co[+2] + XH2(s)': 1,
        'Co[+2] + X(s)': 2,
        'Co[+3] + X(s)': 3,

        'Co(s)': 0,
        'Co[+2]': 1,
        'Co[+3]': 2,
        'Co(HO)2(s)': 3,
        'Co(HO)3(s)': 4,
        'CoHO2[-1]': 5,
        'CoO2(s)': 6,
    }

    sac_mapping = {
        'XCo(s)': 0,
        'XCo(HO)2(s)': 1,
        'XCoO2(s)': 2,
    }
    
    # for entry in vac_entries:
    #     print(entry.name)
    # for entry in sac_entries:
    #     print(entry.name)
    
    for i, entry in enumerate(vac_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        color = vac_colors[vac_mapping[entry.name]]
        ax.fill(x, y, color=color)
    
    for i, entry in enumerate(sac_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        color = sac_colors[sac_mapping[entry.name]]
        ax.fill(x, y, color=color)   
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)

    plt.tight_layout()
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.show()

def main():
    print('\n################## ReCorence Entries ##########################################\n')
    ref_entries = get_ref_entries()
    for entry in ref_entries:
        print(entry)
    
    print('\n################## SAC Entries ##########################################\n')
    sac_entries = get_sac_entries()
    for entry in sac_entries:
        print(entry)

    print('\n################## Solid Entries ##########################################\n')
    solid_entries = get_solid_entries()
    for entry in solid_entries:
        print(entry)
        
    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))
    
    all_entries = ref_entries + sac_entries
    plot_pourbaix(all_entries, f'{filename}_pourbaix_sac{tag}.png')
    
    all_entries = solid_entries + ion_entries
    plot_pourbaix(all_entries, f'{filename}_pourbaix_exp{tag}.png')
    
    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    plot_pourbaix(all_entries, f'{filename}_pourbaix_bulk{tag}.png')


if __name__ == "__main__":
    main()