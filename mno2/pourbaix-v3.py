#!/usr/bin/env python3
"""
Pourbaix Diagram Generation Tool using pymatgen Plotter

This script generates Pourbaix diagrams using pymatgen's PourbaixPlotter,
based on the general structure from pourbaix.py.

Key Features:
- Surface Pourbaix diagram generation from JSON structure files
- Automatic element detection
- Thermodynamic data integration from JSON files
- pymatgen PourbaixPlotter for visualization

Author: Based on pourbaix.py by Hyeonjung Jung
"""

# Standard library imports
import os
import re
import argparse
import json
import glob
import time
from math import log10
from collections import Counter

# Scientific computing imports
import numpy as np
import pandas as pd

# Visualization imports
import matplotlib.pyplot as plt

# Materials science and chemistry imports
from ase.io import read
from mendeleev import element
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

mpr = MPRester('GKTSHzQ7vI80w3OGAA0WThmeiG9OuX5w')
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

def load_jsonc(file_path):
    """Load JSON file with comments (JSONC format)"""
    if not os.path.exists(file_path):
        return {}
    with open(file_path, 'r') as f:
        content = f.read()
    # Remove single-line comments (// ...)
    content = re.sub(r'//.*', '', content)
    return json.loads(content)

# Command line argument parser
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram using pymatgen plotter')
parser.add_argument('--json-dir', type=str, default='.', 
                    help='Folder path containing JSON structure files (default: current folder)')
parser.add_argument('--csv-dir', type=str, default='.', 
                    help='Folder path containing label.csv with species information (default: current folder)')
parser.add_argument('--suffix', type=str, default='', 
                    help='Suffix for the output filename')
parser.add_argument('--label-csv', type=str, 
                    help='Custom path to label.csv file (default: ./label.csv)')
parser.add_argument('--thermo-data', type=str, 
                    help='Custom path to thermodynamic data file (default: ./thermodynamic_data.json)')
parser.add_argument('--ref-energies', type=str, 
                    help='Custom path to reference energies file (default: ./reference_energies.json)')
parser.add_argument('--pHmin', type=float, default=0, 
                    help='Minimum pH value for the diagram (default: 0)')
parser.add_argument('--pHmax', type=float, default=14, 
                    help='Maximum pH value for the diagram (default: 14)')
parser.add_argument('--Umin', type=float, default=-1.5, 
                    help='Minimum potential value in V vs SHE (default: -1.5)')
parser.add_argument('--Umax', type=float, default=+2.0, 
                    help='Maximum potential value in V vs SHE (default: +2.0)')
parser.add_argument('--show-fig', action='store_true', 
                    help='Display the generated plot in a window')
parser.add_argument('--conc', action='store_true', 
                    help='Use concentration for ions')
parser.add_argument('--label', action='store_true', 
                    help='Label domains in the Pourbaix diagram')
parser.add_argument('--HER', action='store_true', 
                    help='Show HER line (0-pH*const)')
parser.add_argument('--OER', action='store_true', 
                    help='Show OER line (1.23-pH*const)')
parser.add_argument('--figx', type=float, default=6, 
                    help='Figure width in inches (default: 6)')
parser.add_argument('--figy', type=float, default=5, 
                    help='Figure height in inches (default: 5)')
parser.add_argument('--exp', action='store_true', 
                    help='Include experimental entries')
parser.add_argument('--simple', action='store_true', 
                    help='Use simple entries')
                    
args = parser.parse_args()

# Physical constants
kjmol = 96.485
calmol = 23.061
kb = 8.617e-5
T = 298.15
const = kb * T * np.log(10)

"""
All computational values from PBE, afm Mn ordering

Error in dH of O2 PBE value:  0.19268103743612963  (if not close to zero: error in O2 PBE energy)
Error in dH of O2 from H2O/H2 PBE values:  0.00019267243612963725  (should be close to zero!)

dH of molecular references at 300K:
1/2 O2 derived from H2O/H2:  -4.658724749999999
1/2 O2 derived from PBE:  -4.851213114999999
1/2 H2 derived from PBE:  -3.2067375399999998
H2O derived from PBE:  -13.578199829999999

Error in dG of O2 PBE value:      0.19014237339508666  (if not close to zero: error in O2 PBE energy)
Error in dG of O2 from H2O/H2 PBE values:  -0.0023459916049133334  (should be close to zero!)

"""

gh2o = -13.578199829999999
gh = -3.2067375399999998
go = -4.658724749999999

TdS_MnO = -0.2314436306161579
TdS_MnO2 = -0.5690729735191998
TdS_MnOOH = -0.8010447219775086
TdS_MnOH2 = -0.8303986596880345
TdS_Mn3O4 = -1.084090519251697
TdS_Mn2O3 = -0.8075002699901537

ZPE_MnO = 0.2867
ZPE_MnO2 = 0.2867
ZPE_MnOOH = 0.5779
ZPE_MnOH2 = 0.8617
ZPE_Mn3O4 = 0.6573
ZPE_Mn2O3 = 0.4937

def main():
    # ========================================
    # SECTION 1: FILE DISCOVERY
    # ========================================
    json_dir = args.json_dir
    json_files = glob.glob(os.path.join(json_dir, "*.json"))
    json_files = sorted(json_files, key=lambda x: os.path.basename(x))
    
    if len(json_files) == 0:
        print("Error: No JSON files found!")
        return
    
    # ========================================
    # SECTION 2: LABEL PROCESSING
    # ========================================
    csv_dir = args.csv_dir
    if hasattr(args, 'label_csv') and args.label_csv:
        label_csv_path = args.label_csv
    else:
        label_csv_path = os.path.join(csv_dir, 'label.csv')
    
    file_labels = {}
    file_oh_counts = {}
    file_gibbs_corrections = {}
    
    if os.path.exists(label_csv_path):
        label_df = pd.read_csv(label_csv_path, header=None, 
                              names=['json_name', 'label', '#OH', 'G_corr'])
        for idx, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']
            if '#OH' in row and not pd.isna(row['#OH']):
                file_oh_counts[json_name] = float(row['#OH'])
            if 'G_corr' in row and not pd.isna(row['G_corr']):
                file_gibbs_corrections[json_name] = float(row['G_corr'])
    else:
        print(f"Warning: Label file {label_csv_path} not found. Using chemical formulas as labels.")
        for json_file in json_files:
            atoms = read(json_file)
            formula = atoms.get_chemical_formula()
            file_labels[os.path.basename(json_file)] = formula
    
    # ========================================
    # SECTION 3: ELEMENT ANALYSIS
    # ========================================
    elements = set()
    for json_file in json_files:
        atoms = read(json_file)
        elements.update(atoms.get_chemical_symbols())
    
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    if 'H' not in sorted_elements:
        sorted_elements.append('H')
    if 'O' not in sorted_elements:
        sorted_elements.append('O')
    sorted_elements = sorted(sorted_elements, key=lambda x: element(x).atomic_number)
    
    # ========================================
    # SECTION 4: BULK DATA CONSTRUCTION FROM JSON FILES
    # ========================================
    print(f"Processing {len(json_files)} bulk structure files...")
    start_time = time.time()
    
    bulk_data = []
    for json_file in json_files:        
        atoms = read(json_file)
        energy = atoms.get_potential_energy()
        symbols = atoms.get_chemical_symbols()
        
        json_basename = os.path.basename(json_file)
        label_name = json_basename.replace('.json', '')
        
        # Apply Gibbs correction if available
        if json_basename in file_gibbs_corrections:
            energy += file_gibbs_corrections[json_basename]
        
        # Get composition from atoms
        comp_dict = Counter(symbols)
        
        # Calculate GCD of all element counts for normalization
        counts = list(comp_dict.values())
        if len(counts) > 0:
            from math import gcd
            from functools import reduce
            gcd_value = reduce(gcd, counts)
        else:
            gcd_value = 1
        
        # Normalize composition and energy by GCD
        normalized_comp_dict = {el: count // gcd_value for el, count in comp_dict.items()}
        normalized_energy = energy / gcd_value
        
        comp_str = ''.join([f"{el}{count}" if count > 1 else el 
                            for el, count in sorted(normalized_comp_dict.items())])
        
        if 'MnO2' in label_name or 'Mn2O4' in label_name or 'Mn4O8' in label_name or 'Mn8O16' in label_name:
            normalized_energy += (ZPE_MnO2 - TdS_MnO2) * normalized_comp_dict['Mn']
        elif 'MnOOH' in label_name:
            normalized_energy += (ZPE_MnOOH - TdS_MnOOH) * normalized_comp_dict['Mn']
        elif 'Mn(OH)2' in label_name:
            normalized_energy += (ZPE_MnOH2 - TdS_MnOH2) * normalized_comp_dict['Mn']
        elif 'MnO' in label_name:
            normalized_energy += (ZPE_MnO - TdS_MnO) * normalized_comp_dict['Mn']
        elif 'Mn3O4' in label_name:
            normalized_energy += (ZPE_Mn3O4 - TdS_Mn3O4) * normalized_comp_dict['Mn'] / 3
        elif 'Mn2O3' in label_name:
            normalized_energy += (ZPE_Mn2O3 - TdS_Mn2O3) * normalized_comp_dict['Mn'] / 2
            
        bulk_data.append({
            'atoms': atoms,
            'energy': normalized_energy,
            'composition': comp_str,
            'name': label_name,
            'json_basename': json_basename
        })
        print(f"composition: {comp_str} (normalized by GCD={gcd_value}), energy: {normalized_energy:.6f} eV/atom, name: {label_name}")
        # print(f"composition: {comp_str}, energy: {energy:.6f} eV/atom, name: {label_name}")
    
    elapsed = time.time() - start_time
    print(f"✓ Structure files processed ({elapsed:.1f} seconds)")
    
    # Determine unique elements for thermodynamic data
    unique_elements = [el for el in sorted_elements if el not in ['H', 'O']]
    
    # ========================================
    # SECTION 8: THERMODYNAMIC DATA LOADING
    # ========================================
    if hasattr(args, 'thermo_data') and args.thermo_data:
        thermo_data_path = args.thermo_data
    else:
        if args.simple:
            thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data_simple.json')
        else:
            thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.json')
        if not os.path.exists(thermo_data_path):
            thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.jsonc')
            
    thermo_data = {}
    if os.path.exists(thermo_data_path):
        thermo_data = load_jsonc(thermo_data_path)
        print(f"Thermodynamic data loaded: {thermo_data_path}")
    
    ref_energies_path = args.ref_energies if hasattr(args, 'ref_energies') and args.ref_energies else os.path.join(os.path.dirname(__file__), 'reference_energies.json')
    ref_energies = {}
    if os.path.exists(ref_energies_path):
        with open(ref_energies_path, 'r') as f:
            ref_energies = json.load(f)
        print(f"Reference energies loaded: {ref_energies_path}")
    
    # ========================================
    # SECTION 9: CONVERT TO PYMATGEN ENTRIES
    # ========================================
    print(f"\n[In progress] Converting to pymatgen PourbaixEntry...")
    
    # Bulk entries from JSON files - convert to formation energy
    bulk_entries = []
    for data in bulk_data:
        try:
            comp_str = data['composition']
            total_energy = data['energy']
            
            # Parse composition to get element counts
            comp = Composition(comp_str)
            
            # Calculate formation energy: E_form = E_total - sum(n_i * E_ref_i)
            # For H and O, use gh and go (chemical potentials from H2O/H2)
            # For other elements, use reference_energies.json
            formation_energy = total_energy
            for el, count in comp.items():
                el_symbol = str(el)
                if el_symbol == 'H':
                    formation_energy -= count * gh
                elif el_symbol == 'O':
                    formation_energy -= count * go
                elif el_symbol in ref_energies:
                    formation_energy -= count * ref_energies[el_symbol]
                else:
                    print(f"Warning: Reference energy for {el_symbol} not found in reference_energies.json")
            
            entry = PourbaixEntry(ComputedEntry(comp_str, formation_energy), entry_id=data['name'], concentration=1)
            bulk_entries.append(entry)
        except Exception as e:
            print(f"Warning: Failed to create entry for {data['name']} (comp: {comp_str}): {e}")

    bulk_entries.append(PourbaixEntry(ComputedEntry('MnOOH', -557.7272/kjmol), entry_id='B-MnOOH_exp', concentration=1))
    bulk_entries.append(PourbaixEntry(ComputedEntry('Mn(OH)2', -615.63376/kjmol), entry_id='D-Mn(OH)2_exp', concentration=1))
    bulk_entries.append(PourbaixEntry(ComputedEntry('Mn3O4', -1283.232/kjmol), entry_id='Mn3O4_exp', concentration=1))
    bulk_entries.append(PourbaixEntry(ComputedEntry('Mn2O3', -881.114/kjmol), entry_id='Mn2O3_exp', concentration=1))
    # bulk_entries.append(PourbaixEntry(ComputedEntry('MnO2', -453.1272/kjmol), entry_id='D-MnO2_exp', concentration=1))

    exp_entries = []
    
    # DFT-calculated data
    # exp_entries.append(PourbaixEntry(ComputedEntry('MnO2', -465.71638408058675/kjmol), entry_id='A-MnO2', concentration=1))
    # exp_entries.append(PourbaixEntry(ComputedEntry('MnOOH', -558.8619486482868/kjmol), entry_id='D-MnOOH', concentration=1))
    # exp_entries.append(PourbaixEntry(ComputedEntry('Mn(OH)2', -604.871812020137/kjmol), entry_id='D-Mn(OH)2', concentration=1))
    # exp_entries.append(PourbaixEntry(ComputedEntry('Mn3O4', -1283.232/kjmol), entry_id='Mn3O4', concentration=1))
    # exp_entries.append(PourbaixEntry(ComputedEntry('Mn2O3', -881.114/kjmol), entry_id='Mn2O3', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('MnO2', -465.138/kjmol), entry_id='B-MnO2_exp', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('MnOOH', -557.7272/kjmol), entry_id='B-MnOOH_exp', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('Mn(OH)2', -615.63376/kjmol), entry_id='D-Mn(OH)2_exp', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('Mn3O4', -1283.232/kjmol), entry_id='Mn3O4_exp', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('Mn2O3', -881.114/kjmol), entry_id='Mn2O3_exp', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('ZnMn8O16', -481.6435639588431*8/kjmol), entry_id='ZnMn8O16', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('ZnMn4O8', -516.324975188399*4/kjmol), entry_id='ZnMn4O8', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('ZnMn2O4', -584.5695268490117*2/kjmol), entry_id='ZnMn2O4', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('ZnMnO2', -584.1993931540368/kjmol), entry_id='ZnMnO2', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('HMn8O16', -472.454106166424*8/kjmol), entry_id='HMn8O16', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('HMn4O8', -482.1357496267614*4/kjmol), entry_id='HMn4O8', concentration=1))
    exp_entries.append(PourbaixEntry(ComputedEntry('HMn2O4', -505.35005103608654*2/kjmol), entry_id='HMn2O4', concentration=1))

    print(f"  Bulk entries: {len(bulk_entries)}")
    
    # Solid entries from thermodynamic data
    solid_entries = []
    for el in unique_elements:
        if el in thermo_data and 'solids' in thermo_data[el]:
            for solid_formula, energy in thermo_data[el]['solids'].items():
                try:
                    energy_ev = energy / calmol
                    entry = PourbaixEntry(PDEntry(solid_formula, energy_ev), concentration=1)
                    solid_entries.append(entry)
                except Exception as e:
                    print(f"Warning: Failed to process solid '{solid_formula}': {e}")

    print(f"  Solid entries: {len(solid_entries)}")
    
    # Ion entries from thermodynamic data
    ion_entries = []
    for el in unique_elements:
        if el in thermo_data and 'ions' in thermo_data[el]:
            for ion_formula, energy in thermo_data[el]['ions'].items():
                try:
                    energy_ev = energy / calmol
                    comp = Ion.from_formula(ion_formula)
                    concentration = 1e-6
                    if args.conc and ion_formula in ['Zn++', 'Mn++']:
                        energy_ev = energy_ev + const * (log10(1) - log10(1e-6))
                        concentration = 1.0
                    entry = PourbaixEntry(IonEntry(comp, energy_ev), concentration=concentration)
                    ion_entries.append(entry)
                except Exception as e:
                    print(f"Warning: Failed to process ion '{ion_formula}': {e}")

    print(f"  Ion entries: {len(ion_entries)}")
    
    # Print all entries list
    print("\n" + "="*70)
    print("All entries list")
    print("="*70)
    
    print(f"\n[Bulk Entries from JSON files] ({len(bulk_entries)}):")
    for i, entry in enumerate(bulk_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV), entry_id: {entry.entry_id}")
    
    print(f"\n[Exp Entries from JSON files] ({len(exp_entries)}):")
    for i, entry in enumerate(exp_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV), entry_id: {entry.entry_id}")
    
    print(f"\n[Solid Entries from thermodynamic_data.json] ({len(solid_entries)}):")
    for i, entry in enumerate(solid_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV)")
    
    print(f"\n[Ion Entries from thermodynamic_data.json] ({len(ion_entries)}):")
    for i, entry in enumerate(ion_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        conc = entry.concentration if hasattr(entry, 'concentration') else 'N/A'
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV, conc: {conc})")
    
    if args.exp:
        all_entries = exp_entries + solid_entries + ion_entries
    else:
        all_entries = bulk_entries + solid_entries + ion_entries
        
    # Remove duplicate entry names, keep only the one with lowest energy
    print(f"\n[In progress] Removing duplicate entry_names...")
    print(f"  Original entry count: {len(all_entries)}")
    
    # Group entries by entry name
    name_groups = {}
    for entry in all_entries:
        entry_name = entry.name if hasattr(entry, 'name') else (entry.entry_id if hasattr(entry, 'entry_id') else 'Unknown')
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        entry_id = entry.entry_id if hasattr(entry, 'entry_id') else 'N/A'
        comp_str = str(entry.composition)

        if entry_name not in name_groups:
            name_groups[entry_name] = []
        name_groups[entry_name].append((entry, energy, entry_id, comp_str))
    
    # Find duplicates and keep only the lowest energy one
    name_dict = {}
    removed_count = 0
    removed_entries = []
    kept_duplicate_entries = []
    
    for entry_name, entries_list in name_groups.items():
        if len(entries_list) > 1:
            # Multiple entries with same name - keep the lowest energy one
            entries_list.sort(key=lambda x: x[1])  # Sort by energy
            
            # Debug output for duplicates
            print(f"\n  Duplicate found: entry_name='{entry_name}' ({len(entries_list)})")
            for i, (entry, energy, entry_id, comp_str) in enumerate(entries_list):
                status = "✓ Keep" if i == 0 else "✗ Remove"
                print(f"    {status}: entry_id={entry_id}, comp={comp_str}, energy={energy:.6f} eV")
            
            kept_entry, kept_energy, kept_entry_id, kept_comp = entries_list[0]
            kept_duplicate_entries.append((entry_name, kept_comp, kept_energy, kept_entry_id))
            
            # Mark removed entries
            for entry, energy, entry_id, comp_str in entries_list[1:]:
                removed_entries.append((entry_name, comp_str, energy, entry_id))
                removed_count += 1
            
            name_dict[entry_name] = kept_entry
        else:
            # Single entry, keep it
            name_dict[entry_name] = entries_list[0][0]
    
    all_entries = list(name_dict.values())
    
    # Print removed entries
    print(f"\n  Removed entry count: {removed_count}")
    print(f"  Remaining entry count: {len(all_entries)}")
    print(f"\n[Total Entries] ({len(all_entries)})")
    print("="*70)
    
    # ========================================
    # SECTION 10: PLOT GENERATION
    # ========================================
    print("\n" + "="*70)
    print("Starting diagram generation")
    print("="*70)
    
    # Combined diagram (all entries)
    print("\n[In progress] Generating Pourbaix diagram...")
    # Add _label suffix if label is enabled
    filename_suffix = args.suffix
    if filename_suffix != '':
        filename_suffix = f'_{filename_suffix}'
    if args.label:
        filename_suffix = f'{filename_suffix}_label' if filename_suffix else '_label'
    plot_pourbaix(all_entries, f'pourbaix{filename_suffix}.pdf', label_domains=args.label, exp_entries=exp_entries, ion_entries=ion_entries)
    
    print("="*70)
    print("All diagrams generated!")
    print("="*70)

def plot_pourbaix(entries, png_name, label_domains=False, exp_entries=None, ion_entries=None):
    """Plot Pourbaix diagram using pymatgen PourbaixPlotter
    
    Args:
        entries: List of PourbaixEntry objects
        png_name: Output filename
        label_domains: If True, label domains in the diagram
        exp_entries: List of experimental entries to group by color
        ion_entries: List of ion entries to extract names from
    """
    if len(entries) == 0:
        print(f"  Warning: Cannot generate diagram - no entries available.")
        return
    
    # Get exp_entry names for color grouping
    exp_entry_names = set()
    if exp_entries:
        for exp_entry in exp_entries:
            exp_entry_names.add(exp_entry.entry_id if hasattr(exp_entry, 'entry_id') and exp_entry.entry_id else exp_entry.name)
    
    # Get ion entry names directly from ion_entries (only Mn species)
    ion_entry_names = set()
    if ion_entries:
        for ion_entry in ion_entries:
            # Get the name from entry (could be entry.name or entry.entry_id)
            ion_name = ion_entry.name if hasattr(ion_entry, 'name') else (ion_entry.entry_id if hasattr(ion_entry, 'entry_id') else str(ion_entry))
            # Only add Mn species to ion_entry_names
            if ion_name.startswith('Mn'):
                ion_entry_names.add(ion_name)
    
    print(f"\n[In progress] Starting {png_name} diagram generation...")
    print(f"  - Entry count: {len(entries)}")
    
    print(f"  - Calculating PourbaixDiagram (convex hull calculation, this may take some time)...")
    start_time = time.time()
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    elapsed = time.time() - start_time
    print(f"  ✓ PourbaixDiagram calculation completed ({elapsed:.1f} seconds)")
    print(f"  - Stable entry count: {len(pourbaix.stable_entries)}")
    
    print(f"  - Initializing Plotter...")
    plotter = PourbaixPlotter(pourbaix)
    
    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    plotter.get_pourbaix_plot(limits=[[args.pHmin, args.pHmax], [args.Umin, args.Umax]], 
                              label_domains=label_domains, label_fontsize=10,
                              show_water_lines=False, show_neutral_axes=False, ax=ax)
    
    stable_entries = pourbaix.stable_entries
    
    # Create mapping from entry names to label.csv names
    # pymatgen formats entry.name with LaTeX subscripts/superscripts (e.g., MnO$_{2}$(s))
    # We need to convert LaTeX format to plain format for matching
    def remove_latex_formatting(text):
        """Remove LaTeX formatting from text"""
        # Remove LaTeX subscripts: $_{...}$
        text = re.sub(r'\$_{([^}]+)}\$', r'\1', text)
        # Remove LaTeX superscripts: $^{...}$ or $^{-...}$
        text = re.sub(r'\$[^{]*\{([^}]+)\}\$', r'\1', text)
        # Remove remaining $ signs
        text = text.replace('$', '')
        return text
    
    # Extract ion patterns from stable entries for color grouping
    # Pattern: Mn[2+], Mn++, Mn[+2], Mn(aq), etc.
    def extract_ion_patterns(entry_name):
        """Extract ion patterns like Mn[2+], Mn++, Zn[2+], etc. from entry name"""
        patterns = []
        plain_name = remove_latex_formatting(entry_name)
        
        # Pattern 1: Mn[2+], Zn[3+], etc. (with brackets)
        matches1 = re.findall(r'([A-Z][a-z]?)\[(\d+)\+?\]', plain_name)
        for element, charge in matches1:
            patterns.append(f"{element}[{charge}+]")
        
        # Pattern 2: Mn++, Zn+++, etc. (multiple + signs without brackets)
        # Match element followed by ++ or +++, etc.
        matches2 = re.findall(r'([A-Z][a-z]?)(\+\+)', plain_name)
        for element, _ in matches2:
            patterns.append(f"{element}[2+]")
        
        # Pattern 3: Mn(aq), Zn(aq), etc.
        matches3 = re.findall(r'([A-Z][a-z]?)\(aq\)', plain_name)
        for element in matches3:
            patterns.append(f"{element}(aq)")
        
        # Pattern 4: Mn$^{2+}$ (LaTeX format in original entry_name)
        matches4 = re.findall(r'([A-Z][a-z]?)\$[^{]*\{(\d+)\+?\}\$', entry_name)
        for element, charge in matches4:
            patterns.append(f"{element}[{charge}+]")
        
        return list(set(patterns))  # Remove duplicates
    
    entry_name_mapping = {}  # Plain format -> entry_id
    entry_latex_mapping = {}  # Plain format -> LaTeX format (for replacement)
    for entry in entries:
        if hasattr(entry, 'entry_id') and entry.entry_id:
            # Convert entry.name from LaTeX to plain format for matching
            plain_name = remove_latex_formatting(entry.name)
            entry_name_mapping[plain_name] = entry.entry_id
            entry_latex_mapping[plain_name] = entry.name  # Store LaTeX format for replacement
    
    # Prepare colors: use Spectral colormap for all entries
    # Group by exp_entry names and ion entry names
    if exp_entry_names or ion_entry_names:
        # Create color map for exp_entries and ion entries using Spectral colormap
        color_groups = {}
        all_groups = list(exp_entry_names) + sorted(list(ion_entry_names))
        n_groups = len(all_groups)
        spectral_colors = [plt.cm.Spectral(i) for i in np.linspace(0.1, 0.9, n_groups)] if n_groups > 0 else []
        for i, group_name in enumerate(all_groups):
            color_groups[group_name] = spectral_colors[i] if spectral_colors else 'silver'

        color_groups = {
            'HMn2O4': 'gold',
            'HMn4O8': 'gold',
            'HMn8O16': 'gold',
            'Mn2HO4': 'gold',
            'Mn4HO8': 'gold',
            'Mn8HO16': 'gold',
            'ZnMnO2': 'greenyellow',
            'ZnMn2O4': 'greenyellow',
            'ZnMn4O8': 'greenyellow',
            'ZnMn8O16': 'greenyellow',
            'MnZnO2': 'greenyellow',
            'Mn2ZnO4': 'greenyellow',
            'Mn4ZnO8': 'greenyellow',
            'Mn8ZnO16': 'greenyellow',
            'Mn[+2]': 'white',
            'Mn[+3]': 'white',
            'MnO4[-1]': 'white',
            'MnO4[-2]': 'white',
            'MnHO2[-1]': 'white',
            'Mn(OH)2': 'lightblue',
            'Mn(HO)2': 'lightblue',
            'MnOOH': 'lightblue',
            'MnO2': 'mediumpurple',
            'Mn2O3': 'silver',
            'Mn3O4': 'silver',
            'Mn': 'silver',
        }        
        # for color_group in color_groups:
        #     print(f"  - {color_group}: {color_groups[color_group]}")

        # Create color mapping: entries containing same Mn species get same color
        # Don't distinguish between exp_entry and ion_entry - prioritize Mn species
        entry_color_map = {}
        unmatched_entries = []
        
        # Sort stable_entries by energy for consistent processing
        stable_entries_sorted = sorted(stable_entries, key=lambda e: e.energy if hasattr(e, 'energy') else (e.entry.energy if hasattr(e, 'entry') and hasattr(e.entry, 'energy') else 0))
        
        # Convert color_groups keys to list for iteration (preserves insertion order)
        # This ensures we check color groups in the order they are defined above
        color_group_list = list(color_groups.keys())

        # Single pass: process all entries with priority (Mn species > exp_entry)
        for entry in stable_entries_sorted:
            entry_name = entry.name            
            matched_group = None
            
            for color_group in color_group_list:
                # print(f"  - {color_group} in {entry_name}")
                if color_group in entry_name:
                    matched_group = color_group
                    break
            # print(f"  - {entry_name} matched to {matched_group}")

            # Assign color from color_groups if matched
            if matched_group and matched_group in color_groups:
                entry_color_map[entry_name] = color_groups[matched_group]
            else:
                # Track unmatched entries to assign colors later
                unmatched_entries.append(entry_name)
        
        # Assign Spectral colors to unmatched entries
        n_unmatched = len(unmatched_entries)
        if n_unmatched > 0:
            unmatched_colors = [plt.cm.Grays(i) for i in np.linspace(0.1, 0.9, n_unmatched)]
            for i, entry_name in enumerate(unmatched_entries):
                # entry_color_map[entry_name] = unmatched_colors[i]
                entry_color_map[entry_name] = 'silver'
    else:
        # Use Spectral colormap for all entries
        n_entries = len(stable_entries)
        spectral_colors = [plt.cm.Spectral(i) for i in np.linspace(0.1, 0.9, n_entries)] if n_entries > 0 else []
        
        # Create color mapping for each entry
        entry_color_map = {}
        for i, entry in enumerate(stable_entries):
            entry_color_map[entry.name] = spectral_colors[i] if spectral_colors else 'silver'
    
    # Fill each domain with color
    for entry in stable_entries:
        try:
            vertices = plotter.domain_vertices(entry)
            if len(vertices) > 0:
                x, y = zip(*vertices)
                color = entry_color_map.get(entry.name, 'silver')
                ax.fill(x, y, color=color, alpha=0.6, edgecolor='black', linewidth=0.5)
        except Exception as e:
            print(f"Warning: Could not fill domain for {entry.name}: {e}")
    
    # Style adjustments and label replacement
    for line in ax.lines:
        line.set_linewidth(0.5)
    for text in ax.texts:
        # Reduce font size
        text.set_fontsize(10)
        text.set_color('black')
        text.set_fontweight('bold')
        
        # Replace entry name with label.csv name if available
        old_text = text.get_text()
        
        # Convert old_text from LaTeX format to plain format for matching
        plain_text = remove_latex_formatting(old_text)
        
        # Try to match with entry names in plain format
        for plain_name, label_name in entry_name_mapping.items():
            if plain_name in plain_text:
                # Get the LaTeX format for this entry
                latex_name = entry_latex_mapping.get(plain_name, plain_name)
                
                # Replace LaTeX formatted name with label_name in old_text
                if latex_name in old_text:
                    new_text = old_text.replace(latex_name, label_name)
                    text.set_text(new_text)
                    break
                elif plain_name == plain_text:
                    # If exact match and LaTeX not found, replace entire text
                    text.set_text(label_name)
                    break
    
    pH2 = np.arange(args.pHmin, args.pHmax + 0.01, 0.01)
    if args.OER:
        plt.plot(pH2, 1.23 - pH2 * const, linestyle='--', color='red', lw=2.0)
    if args.HER:
        plt.plot(pH2, 0 - pH2 * const, linestyle='--', color='red', lw=2.0)
    plt.plot(pH2, 2.6 -0.76 -1*const -pH2 * const, linestyle='--', color='blue')
    plt.plot(pH2, 2.3 -0.76 -1*const -pH2 * const, linestyle='--', color='blue')

    
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    print(f"  - Saving plot...")
    plt.tight_layout()
    plt.savefig(png_name, dpi=200, bbox_inches='tight', transparent=True)
    if args.show_fig:
        plt.show()
    plt.close()
    
    total_time = time.time() - start_time
    print(f"  ✓ {png_name} diagram generation completed! (Total {total_time:.1f} seconds)\n")

if __name__ == "__main__":
    main()
