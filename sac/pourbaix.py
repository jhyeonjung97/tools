#!/usr/bin/env python3

import os
import sys
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.entries.computed_entries import ComputedEntry

# Ignore warning messages
warnings.filterwarnings('ignore')

# Set plot file name
PLOT_NAME = 'Pourbaix_PBE+U_U2p75_allLayered'

# Load Materials Project API key from environment variable
API_KEY = os.getenv('MAPI_KEY')
if not API_KEY:
    sys.exit("Error: MAPI_KEY environment variable not set.")

mpr = MPRester(API_KEY)

# Define constants
kJmol = 96.485

# gas
h2 = -6.77149190
zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408
gh2 = h2 + zpeh2 - tsh2 + cvh2

elements_data = {
    "Ti": {"electrode_potential": -1.63, "cation_charge": 2},  # Ti^2+ + 2e- → Ti (Ti^2+ prioritized over Ti^4+)
    "V":  {"electrode_potential": -1.18, "cation_charge": 2},  # V^2+ + 2e- → V (V^2+ prioritized over V^3+)
    "Cr": {"electrode_potential": -0.91, "cation_charge": 2},  # Cr^2+ + 2e- → Cr (Cr^2+ prioritized over Cr^3+)
    "Mn": {"electrode_potential": -1.18, "cation_charge": 2},  # Mn^2+ + 2e- → Mn
    "Fe": {"electrode_potential": -0.44, "cation_charge": 2},  # Fe^2+ + 2e- → Fe (Fe^2+ prioritized over Fe^3+)
    "Co": {"electrode_potential": -0.28, "cation_charge": 2},  # Co^2+ + 2e- → Co (Co^2+ prioritized over Co^3+)
    "Ni": {"electrode_potential": -0.25, "cation_charge": 2},  # Ni^2+ + 2e- → Ni
    "Cu": {"electrode_potential": 0.34,  "cation_charge": 2},  # Cu^2+ + 2e- → Cu

    "Zr": {"electrode_potential": -1.45, "cation_charge": 2},  # Zr^2+ + 2e- → Zr (Zr^2+ prioritized over Zr^4+)
    "Nb": {"electrode_potential": -1.00, "cation_charge": 3},  # Nb^3+ + 3e- → Nb (Nb^3+ prioritized over Nb^5+)
    "Mo": {"electrode_potential": -0.20, "cation_charge": 3},  # Mo^3+ + 3e- → Mo (Mo^3+ prioritized over Mo^6+)
    "Tc": {"electrode_potential": 0.74,  "cation_charge": 7},  # Tc^7+ + 7e- → Tc (No lower oxidation state available)
    "Ru": {"electrode_potential": 0.45,  "cation_charge": 2},  # Ru^2+ + 2e- → Ru (Ru^2+ prioritized over Ru^3+)
    "Rh": {"electrode_potential": 0.76,  "cation_charge": 2},  # Rh^2+ + 2e- → Rh (Rh^2+ prioritized over Rh^3+)
    "Pd": {"electrode_potential": 0.95,  "cation_charge": 2},  # Pd^2+ + 2e- → Pd

    "Hf": {"electrode_potential": -1.55, "cation_charge": 4},  # Hf^4+ + 4e- → Hf (No lower oxidation state available)
    "Ta": {"electrode_potential": -1.10, "cation_charge": 3},  # Ta^3+ + 3e- → Ta (Ta^3+ prioritized over Ta^5+)
    "W":  {"electrode_potential": -0.11, "cation_charge": 4},  # W^4+ + 4e- → W (W^4+ prioritized over W^6+)
    "Re": {"electrode_potential": 0.30,  "cation_charge": 4},  # Re^4+ + 4e- → Re (Re^4+ prioritized over Re^7+)
    "Os": {"electrode_potential": 0.40,  "cation_charge": 2},  # Os^2+ + 2e- → Os (Os^2+ prioritized over Os^4+)
    "Ir": {"electrode_potential": 1.16,  "cation_charge": 2},  # Ir^2+ + 2e- → Ir (Ir^2+ prioritized over Ir^3+)
    "Pt": {"electrode_potential": 1.20,  "cation_charge": 2},  # Pt^2+ + 2e- → Pt
}

potential = elements_data['Fe']['electrode_potential']
charge = elements_data['Fe']['cation_charge']

# metal_path = '/pscratch/sd/j/jiuy97/6_MNC/gas/metals.tsv'
metal_path = '/Users/hailey/Desktop/jan20/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
metal = metal_df.loc['Fe', 'energy']
N4C26 = -.27195317E+03

# Read the TSV file with correct delimiter and index column
# df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/1Fe_energies.tsv', delimiter='\t', index_col=0)
df = pd.read_csv('/Users/hailey/Desktop/jan20/1Fe_energies.tsv', delimiter='\t', index_col=0)

# Process the composition column
df['name'] = 'FeNC(' + df.index.str.upper() + ')'
df['comp'] = 'FeHH' + df.index.str.upper().str.replace("-", "")
df['comp'] = df['comp'].str.replace('FeHHVAC', 'VAC')
df['comp'] = df['comp'].str.replace('FeHHCLEAN', 'FeHH')
df['comp'] = df['comp'].str.replace('FeHHMH', 'FeHHH')
df['comp'] = df['comp'].str.replace('FeHHNH', 'FeHHH')

# Assuming charge, potential, and gh2 variables are defined elsewhere
df['energy'] = df['dG'] + df.loc['clean', 'G'] - metal - N4C26
df = df.drop(index='vac')
df = df.dropna()

def get_solid_entries():
    """Generate solid entries."""    
    solid_entries = []
    
    for index, row in df.iterrows():
        try:
            comp = Ion.from_formula(row['comp'])
            energy = row['energy']
            print(energy)
            name = row['name']
            
            # Create IonEntry with zero charge by default
            entry = PourbaixEntry(IonEntry(comp, energy, name=name))
            solid_entries.append(entry)
        
        except Exception as e:
            print(f"Error processing entry {row['comp']}: {e}")
            
    return solid_entries

def get_ion_entries():
    """Fetch ion entries from the Materials Project API."""
    try:
        ion_entries = mpr.get_pourbaix_entries(["Fe"])
    except Exception as e:
        sys.exit(f"Error fetching ion data from Materials Project API: {e}")

    return ion_entries

def plot_pourbaix(entries):
    """Plot and save Pourbaix diagram."""
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    # Generate the Pourbaix plot and get the Axes object
    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-2, 4]])
    
    # Customize the plot
    for line in ax.lines:
        line.set_linewidth(1.0)  # Adjust line thickness
    
    for text in ax.texts:
        text.set_fontsize(14)  # Adjust phase label font size
    
    # Set axis labels and tick font size
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    fig = ax.figure
    fig.set_size_inches((8, 7))
    plt.tight_layout()

    plt.savefig(PLOT_NAME + '.pdf')
    plt.savefig(PLOT_NAME + '.png')
    plt.show()

def main():
    """Main execution function."""
    print('\n################## Solid Entries ##########################################\n')
    solid_entries = get_solid_entries()
    for entry in solid_entries:
        print(entry)

    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = solid_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))
    
    plot_pourbaix(all_entries)

if __name__ == "__main__":
    main()