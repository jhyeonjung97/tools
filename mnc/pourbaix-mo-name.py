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

png_name = '3Mo_pourbaix'
tsv_name = '3Mo_energies.tsv'
        
kJmol = 96.485
calmol = 23.061
water = 2.4583 # the standard Gibbs free energy of formation of water

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.560
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

metal_path = './metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
gm = metal_df.loc['Mo', 'energy']

df = pd.read_csv(tsv_name, delimiter='\t', index_col=0)

df['name'] = 'MoNC(' + df.index.str.upper() + ')'
df['comp'] = 'MoX' + df.index.str.upper().str.replace("-", "")
df['comp'] = df['comp'].str.replace('MoXVAC', 'H2X')
df['name'] = df['name'].str.replace('MoNC(VAC)', 'Mo⁺²+H₂NC', regex=False)
df['comp'] = df['comp'].str.replace('MoXCLEAN', 'MoX')
df['name'] = df['name'].str.replace('MoNC(CLEAN)', 'MoNC(clean)')
df['comp'] = df['comp'].str.replace('MoXMH', 'MoXH')
df['comp'] = df['comp'].str.replace('MoXNH', 'MoXH')

# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - 2 * gn2 - 26 * gc - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - N4C26 - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
df['energy'] = df['dG'] + df.loc['clean', 'G'] + gh2 - gm - H2N4C26 - 2 * dgh - water * (df['#O'] + df['#OH'] + df['#OOH']*2)

df = df.drop(index='vac')
df = df.drop(index='o-oh')
df = df.drop(index='ooho')
df = df.drop(index='o-ooh')
df = df.drop(index='ooh-o')
df = df.drop(index='oh-ooh')
df = df.drop(index='ooh-oh')
df = df.drop(index='ooh-ooh')
df = df.dropna()
print(df)

def get_ref_entries():
    ref_entries = []
    
    refs={
        'Mo': 'Mo(s)',
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
        'Mo': 0,
        'MoO2' : -120.000/calmol,
        # 'MoO3': -227.000/calmol,
        'MoO3': -161.950/calmol,
        }
    
    for solid, energy in solids.items():
        entry = PourbaixEntry(PDEntry(solid, energy))
        solid_entries.append(entry)

    return solid_entries

def get_ion_entries():
    ion_entries = []
    ions={
        'Mo+++': -13.800/calmol,
        # 'HMoO4-': -213.600/calmol,
        'MoO4--': -205.420/calmol,
        }
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy), concentration=1e-6)
        ion_entries.append(entry)

    return ion_entries
    
def plot_pourbaix(entries, png_name):    
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    fig, ax = plt.subplots(figsize=(6, 5))    
    plotter.get_pourbaix_plot(limits=[[0, 14], [-1, 3]], label_domains=False, label_fontsize=14,
                              show_water_lines=False, show_neutral_axes=False, ax=ax)
    stable_entries = pourbaix.stable_entries

    for line in ax.lines:
        line.set_linewidth(0.5)
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_color('black')
        text.set_fontweight('bold')
    
    if 'sac' in png_name:
        ax.text(0.2, -0.9, r"S$_{\mathbf{0}}$+Mo(s)", fontsize=14, color="black", fontweight='bold')
        ax.text(7.0, 1.5, r"S$_{\mathbf{9}}$", fontsize=14, color="black", fontweight='bold', ha='center', va='center')
        ax.text(7.0, -0.25, r"S$_{\mathbf{4}}$", fontsize=14, color="black", fontweight='bold', ha='center', va='center')
    elif 'bulk' in png_name:
        ax.text(0.2, -0.9, r"S$_{\mathbf{0}}$+Mo(s)", fontsize=14, color="black", fontweight='bold')
        ax.text(3.0, -0.3, r"S$_{\mathbf{0}}$+MoO$_{\mathbf{2}}$(s)", fontsize=14, color="black", fontweight='bold')
        ax.text(0.2, 0.0, r"S$_{\mathbf{0}}$+Mo$^{\mathbf{3+}}$", fontsize=14, color="black", fontweight='bold')
        ax.text(8.0, -0.1, r"S$_{\mathbf{0}}$+MoO$_{\mathbf{4}}^{\mathbf{2-}}$", fontsize=14, color="black", fontweight='bold')
        ax.text(0.2, 0.7, r"S$_{\mathbf{0}}$+MoO$_{\mathbf{3}}$(s)", fontsize=14, color="black", fontweight='bold')
        ax.text(8.0, 1.9, r"S$_{\mathbf{v}}$+MoO$_{\mathbf{4}}^{\mathbf{2-}}$", fontsize=14, color="black", fontweight='bold')
        ax.text(0.2, 1.9, r"S$_{\mathbf{v}}$+MoO$_{\mathbf{3}}$(s)", fontsize=14, color="black", fontweight='bold')

    vac_entries = [entry for entry in stable_entries if 'XMo' not in entry.name]
    sac_entries = [entry for entry in stable_entries if 'XMo' in entry.name]
    vac_colors = [plt.cm.Greys(i) for i in np.linspace(0.1, 0.3, len(vac_entries))]
    sac_colors = [plt.cm.Purples(i) for i in np.linspace(0.1, 0.3, len(sac_entries))]
    
    vac_mapping = {
        'Mo(s) + XH2(s)': 0,
        'MoO2(s) + XH2(s)': 1,
        'Mo[+3] + XH2(s)': 2,
        'MoO4[-2] + XH2(s)': 3,
        'MoO3(s) + XH2(s)': 4,
        'MoO4[-2] + X(s)': 5,
        'MoO3(s) + X(s)': 6,
        'XH2(s) + Mo(s)': 0,
        'XH2(s) + MoO2(s)': 1,
        'XH2(s) + Mo[+3]': 2,
        'XH2(s) + MoO4[-2]': 3,
        'XH2(s) + MoO3(s)': 4,
        'X(s) + MoO4[-2]': 5,
        'X(s) + MoO3(s)': 6,
    }
    
    sac_mapping = {
        'XMoO(s)': 0,
        'XMoO2(s)': 1,
    }
    
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
    plt.savefig(png_name, dpi=300, bbox_inches='tight')
    plt.show()

def main():
    print('\n################## ReMorence Entries ##########################################\n')
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
    plot_pourbaix(all_entries, f'{png_name}_sac_name.png')
    
    # plot_pourbaix(solid_entries, f'{png_name}_solid.png')
    # plot_pourbaix(ion_entries, f'{png_name}_ion.png')
    # all_entries = solid_entries + ion_entries
    # plot_pourbaix(all_entries, f'{png_name}_exp.png')

    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    plot_pourbaix(all_entries, f'{png_name}_bulk_name.png')


if __name__ == "__main__":
    main()