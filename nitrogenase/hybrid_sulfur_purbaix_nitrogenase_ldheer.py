#!/usr/bin/env python3

import os
import sys
import warnings
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Set CMU Sans Serif font globally if available.
# If the font file is not present, matplotlib will fall back to its default font.
font_path = '/Users/ldheer/Library/Fonts/cmunss.otf'
if os.path.exists(font_path):
    cmu_sans = fm.FontProperties(fname=font_path)
    font_family = cmu_sans.get_name()
else:
    font_family = 'DejaVu Sans'

plt.rcParams.update({
    'font.family': font_family,
    #'font.sans-serif': ['CMU Sans Serif'],  # You must have this font installed
    'pdf.fonttype': 42,  # Ensure text is editable in PDF
    'ps.fonttype': 42
})

from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

warnings.filterwarnings('ignore')
png_name = 'nitrogenase'
tag = 'BS7_U2_DFT_vac_new'

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
tsh2s_aq = 0.93 # RPBE, aq 10-6 M (H2S)
tsh2s_gas = 0.636 # RPBE, aq

dgh2s_aq = zpeh2s - tsh2s_aq + cvh2s
dgh2s_gas = zpeh2s - tsh2s_gas + cvh2s
gh2s_aq = h2s + zpeh2s - tsh2s_aq + cvh2s
gh2s_gas = h2s + zpeh2s - tsh2s_gas + cvh2s

# ------------------------------------------------------------------
# DFT sulfur reference used ONLY for the FeMo-cofactor E-state ladder.
# This keeps E0/E1/E2/E3/E4 relative to V_S consistent with the
# original DFT thermodynamic scheme.
# ------------------------------------------------------------------
s = -126.00103840/32 # RPBE; not used below, printed only for reference
s_h2s = h2s - h2
print(s, s_h2s, s_h2s+0.183)
# Use the original variable names downstream for minimal code changes.
gs = s_h2s+0.183-0.1
#Formation energy for H2S with reference to S solid (0.183)
#Entropy difference between solid S and H2S S (-0.1)
fh2s_aq = gh2s_aq - gh2 - gs
fh2s_gas = gh2s_gas - gh2 - gs
print('dgh2s_aq=', dgh2s_aq)
print('s_32=', s)
print('DFT sulfur reference from H2S-H2 =', gs)



# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

dgoh = zpeoh + cvoh - tsoh
dgo = zpeo + cvo - tso
dgh = 0.25
dgs = 0.05


df = pd.DataFrame()
df.loc['EO', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-166.574600, 2, 8, 7, 1, 1, 2, 10]
df.loc['E1', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-170.406509, 2, 9, 7, 1, 1, 2, 10]
df.loc['E2', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-174.416577, 2, 10, 7, 1, 1, 2, 10]
df.loc['E3', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-178.024417, 2, 11, 7, 1, 1, 2, 10]
#df.loc['E4', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-181.091342+0.14, 2, 12, 7, 1, 1, 2, 10]
df.loc['Vac', ['E', '#C', '#H', '#Fe', '#Mo', '#N', '#O', '#S']] = [-168.645900+0.14, 2, 10, 7, 1, 1, 2, 9] # C2H8Fe7MoNO2S9 + H2

df['comp'] = (
    'X'
    + 'H' + (df['#H']-8).astype(str)
    + 'S' + (df['#S']-9).astype(str)
)
df['G'] = df['E'] + (0 + dgh * (df['#H']-8)) + (0 + dgs * (df['#S']-9)) - (0 + gh * (df['#H']-8)  + gs * (df['#S']-9))
df['energy'] = df['G'] - df.loc['Vac', 'G'] # - water * df['#O']
print(df)

def get_ref_entries():
    """No extra reference entries are needed for this trial hybrid diagram.

    The FeMo-cofactor states are included as cluster entries. The dissolved
    sulfur species H2S(aq), HS-, and S2- are included as aqueous ion entries.
    """
    return []
    
def get_cluster_entries():
    """DFT-derived E-state and V_S cluster entries."""
    cluster_entries = []
    
    for index, row in df.iterrows():    
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy']))
        cluster_entries.append(entry)
    
    return cluster_entries  
          
def get_gas_entries():
    """No gas-phase sulfur species in this trial.

    This keeps the vacancy-associated sulfur chemistry aqueous only:
    V_S + H2S(aq), V_S + HS-, and V_S + S2-.
    """
    return []
          
def get_ion_entries(concentration=1e-6):
    """Experimental aqueous sulfide entries.

    These values are used only to generate the pH-dependent dissolved sulfur
    variants of the vacancy state:
        V_S + H2S(aq)
        V_S + HS-
        V_S + S2-

    Values are in kcal/mol and converted to eV using:
        1 eV = 23.061 kcal/mol
    """
    ion_entries = []
    kcalmol = 23.061

    ions = {
        #'H2S': -7.892 / kcalmol,  # H2S(g)
        'HS-': 3.010 / kcalmol,
        'S--': 21.958 / kcalmol,
        'H2S': -6.540 / kcalmol,  #H2S(aq)
        'HS-': 3.010 / kcalmol,
        'S--': 21.958 / kcalmol,
        'S2--': 19.749 / kcalmol,
        'S3--': 17.968 / kcalmol,
        'S4--': 16.615 / kcalmol,
        'S5--': 15.689 / kcalmol,
        'H2S2O3': -129.900 / kcalmol,
        'HS2O3-': -129.500 / kcalmol,
        'S2O3--': -127.200 / kcalmol,
        'S5O6--': -228.500 / kcalmol,
        'S4O6--': -244.300 / kcalmol,
        'HS2O4-': -141.408 / kcalmol,
        'S2O4--': -138.000 / kcalmol,
        'S3O6--': -229.000 / kcalmol,
        'H2SO3': -128.690 / kcalmol,
        'HSO3-': -126.000 / kcalmol,
        'SO3--': -116.100 / kcalmol,
        'S2O6--': -231.000 / kcalmol,
        'H2SO4': -177.340 / kcalmol,
        'HSO4-': -179.940 / kcalmol,
        'SO4--': -177.340 / kcalmol,
        'S2O8--': -262.000 / kcalmol,
    }
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy))
        print(ion, energy/comp['S'])
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
        
    pH2 = np.arange(0, 14.01, 0.01)
    plt.plot(pH2, 0.05 - pH2 * const, color='#55C5FF', lw=3, dashes=(4, 3))
    plt.plot(pH2, -0.7 - pH2 * const, color='#E30202', lw=3, dashes=(4, 3))
    ax.text(2.5, -1.2, r"$\mathrm{E^{0}_{\left( ATP / ADP \right)}}$ = -0.7 V vs RHE", 
            color='#E30202', rotation=-14, fontsize=16)
    ax.text(2.5, -0.45, r"$\mathrm{E^{0}_{\left( N_{2} / NH_{3} \right)}}$ = 0.05 V vs RHE",
            color='#55C5FF', rotation=-14, fontsize=16)
   
    name_mapping1 = {
        'H2S(aq) + XH2(s)': 'XH2(s) + H2S(aq)',
        'H2S(s) + XH2(s)': 'XH2(s) + H2S(aq)',
        'HS[-1] + XH2(s)': 'XH2(s) + HS[-1]',
        'S[-2] + XH2(s)': 'XH2(s) + S[-2]',
        'H2SO4(aq) + XH2(s)': 'XH2(s) + H2SO4(aq)',
        'HSO4[-1] + XH2(s)': 'XH2(s) + HSO4[-1]',
        'SO4[-1] + XH2(s)': 'XH2(s) + SO4[-1]',
        'SO4[-2] + XH2(s)': 'XH2(s) + SO4[-2]',
    }
    
    name_mapping2 = {
        'XS(s)': r'E$_{0}$', 
        'XHS(s)': r'E$_{1}$',
        'XH2S(s)': r'E$_{2}$',
        'XH3S(s)': r'E$_{3}$',
        'XH4S(s)': r'E$_{4}$',
        'XH2(s)': r'E$_{v}$',
        'S[-2]': 'S²⁻',
        'HS[-1]': 'HS⁻',
        'H2S(s)': 'H₂S(aq)',
        'H2S(aq)': 'H₂S(aq)',
        'SO4[-1]': ' S₂O₈²⁻',
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
        'H2S(aq)',
        'H2S(s)',
        'HS[-1]',
        'S[-2]',
    ]
        
    sulfur_entries = [entry for entry in stable_entries if 'XH2(s)' not in entry.name]
    vacancy_entries = [entry for entry in stable_entries if 'XH2(s)' in entry.name]
    sulfur_colors = {
    'XS(s)': (1.0, 1.0, 1.0),  # E0: White
    'XHS(s)': (230/255, 230/255, 230/255),  # E2: Light Grey
    'XH2S(s)': (220/255, 220/255, 220/255),  # E3: Grey
    'XH3S(s)': (175/255, 175/255, 175/255),  # Additional sulfur state
    'XH4S(s)': (175/255, 150/255, 175/255)   # Additional sulfur state
}
    vacancy_colors = {
    'H2S(aq)':   (181/255, 138/255, 88/255),
    'H2S(s)':    (181/255, 138/255, 88/255),
    'HS[-1]':    (203/255, 173/255, 138/255),
    'S[-2]':     (234/255, 222/255, 208/255),
    'HSO4[-1]':  (250/255, 172/255, 212/255),
    'SO4[-1]':   '#C8D8BF',
    'H2SO4(aq)': '#C8D8BF',
    'SO4[-2]':   (252/255, 213/255, 212/255),

} # sulfide-only vacancy variants
        
    for i, entry in enumerate(sulfur_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        color = sulfur_colors.get(entry.name, '#FFFFFF')
        ax.fill(x, y, color=color)
    
    for i, entry in enumerate(vacancy_entries):
        vertices = plotter.domain_vertices(entry)
        x, y = zip(*vertices)
        ion_part = entry.name.replace(' + ', '').replace('XH2(s)', '')
        color = vacancy_colors.get(ion_part, '#FFFFFF')
        ax.fill(x, y, color=color)
        
    if 'bulk' in png_name:
        ax.text(12, -0.50, r"$\mathrm{E}_{0}$", fontsize=18, ha='center', fontweight='bold')
        ax.text(12, -0.75, r"$\mathrm{E}_{2}$", fontsize=18, ha='center', fontweight='bold')
        ax.text(12, -1.00, r"$\mathrm{E}_{3}$", fontsize=18, ha='center', fontweight='bold')
        ax.text(3.5, -1.6, r"$\mathrm{E}_{4}\mathrm{(*4HS)}$", fontsize=18, ha='center', fontweight='bold')
        ax.text(2.8, -1.6, r"$\mathrm{V_S} + \mathrm{H_2S(aq)}$", fontsize=18, ha='center', fontweight='bold')
        ax.text(9.0, -1.6, r"$\mathrm{V_S} + \mathrm{HS}^{-}$", fontsize=18, fontweight='bold')
        ax.text(13.8, -1.9, r"$\mathrm{V_S} + \mathrm{S}^{2-}$", fontsize=18, ha='right', fontweight='bold')
        ax.text(0.7, 0.8, r"$\mathrm{V_S} + \mathrm{HSO}_{4}^{-}$", fontsize=18, fontweight='bold')
        ax.text(9.0, 0.5, r"$\mathrm{V_S} + \mathrm{SO}_{4}^{2-}$", fontsize=18, fontweight='bold')

        #ax.text(12, -0.50, r"$E_{0}$", fontsize=14, ha='center')"
        #ax.text(12, -0.75, r"$E_{2}$", fontsize=14, ha='center')
        #ax.text(12, -1.00, r"$E_{3}$", fontsize=14, ha='center')
        #ax.text(0.7, 0.8, r"$E_{v} + HSO_{4}^{-}$", fontsize=14)
        #ax.text(9.0, 0.5, r"$E_{v} + SO_{4}^{2-}$", fontsize=14)
        #ax.text(3.5, -1.6, r"$E_{v} + H_{2}S(aq)$", fontsize=14, ha='center')
        #ax.text(9.0, -1.6, r"$E_{v} + HS^{-}$", fontsize=14)
        #ax.text(13.8, -1.9, r"$E_{v} + S^{2-}$", fontsize=14, ha='right')

        #ax.text(12, -0.50, r"$\mathbf{E_{0}}$", fontsize=18, ha='center')
        #ax.text(12, -0.75, r"$\mathbf{E_{2}}$", fontsize=18, ha='center')
        #ax.text(12, -1.00, r"$\mathbf{E_{3}}$", fontsize=18, ha='center')
        #ax.text(0.7, 0.8, r"$\mathbf{V_S} + \mathbf{HSO}_{4}^{-}$", fontsize=18)
        #ax.text(9.0, 0.5, r"$\mathbf{V_S} + \mathbf{SO}_{4}^{2-}$", fontsize=18)
        #ax.text(3.5, -1.6, r"$\mathbf{V_S} + \mathbf{H_{2}S(aq)}$", fontsize=18, ha='center')
        #ax.text(9.0, -1.6, r"$\mathbf{V_S} + \mathbf{HS}^{-}$", fontsize=18)
        #ax.text(13.8, -1.9, r"$\mathbf{V_S} + \mathbf{S}^{2-}$", fontsize=18, ha='right') 
        
    ax.set_xlabel("pH", fontsize=16)
    ax.set_ylabel("Potential (V)", fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    ax.text(0.99, 0.99, r"$\mathrm{GGA}+U_{\mathrm{Fe}}\;(2 \, \mathrm{eV})$", 
        transform=ax.transAxes,  # Position relative to axes
        fontsize=16, 
        fontweight='bold', 
        ha='right',  # Align text to the right
        va='top')    # Align text to the top    
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
    plot_pourbaix(all_entries, f'{png_name}_bulk{tag}.pdf')
    
if __name__ == "__main__":
    main()
