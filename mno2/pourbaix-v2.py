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
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

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
parser.add_argument('--Umin', type=float, default=-1, 
                    help='Minimum potential value in V vs SHE (default: -1)')
parser.add_argument('--Umax', type=float, default=3, 
                    help='Maximum potential value in V vs SHE (default: 3)')
parser.add_argument('--show-fig', action='store_true', 
                    help='Display the generated plot in a window')
parser.add_argument('--conc', action='store_true', 
                    help='Use concentration for ions')
                    
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

MnO_TdS=  -0.2314436306161579
MnO2_TdS=  -0.5690729735191998
MnOOH_TdS=  -0.8010447219775086
MnOH2_TdS=  -0.8303986596880345
Mn3O4_TdS=  -1.084090519251697
Mn2O3_TdS=  -0.8075002699901537

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
            normalized_energy -= MnO2_TdS * normalized_comp_dict['Mn']
        elif 'MnOOH' in label_name:
            normalized_energy -= MnOOH_TdS * normalized_comp_dict['Mn']
        elif 'Mn(OH)2' in label_name:
            normalized_energy -= MnOH2_TdS * normalized_comp_dict['Mn']
        elif 'MnO' in label_name:
            normalized_energy -= MnO_TdS * normalized_comp_dict['Mn']
        elif 'Mn3O4' in label_name:
            normalized_energy -= Mn3O4_TdS * normalized_comp_dict['Mn'] / 3
        elif 'Mn2O3' in label_name:
            normalized_energy -= Mn2O3_TdS * normalized_comp_dict['Mn'] / 2
            
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
    print(f"✓ 구조 파일 처리 완료 ({elapsed:.1f}초 소요)")
    
    # Determine unique elements for thermodynamic data
    unique_elements = [el for el in sorted_elements if el not in ['H', 'O']]
    
    # ========================================
    # SECTION 8: THERMODYNAMIC DATA LOADING
    # ========================================
    if hasattr(args, 'thermo_data') and args.thermo_data:
        thermo_data_path = args.thermo_data
    else:
        thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.json')
        if not os.path.exists(thermo_data_path):
            thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.jsonc')
    
    thermo_data = {}
    if os.path.exists(thermo_data_path):
        thermo_data = load_jsonc(thermo_data_path)
        print(f"Thermodynamic data 로드 완료: {thermo_data_path}")
    
    ref_energies_path = args.ref_energies if hasattr(args, 'ref_energies') and args.ref_energies else os.path.join(os.path.dirname(__file__), 'reference_energies.json')
    ref_energies = {}
    if os.path.exists(ref_energies_path):
        with open(ref_energies_path, 'r') as f:
            ref_energies = json.load(f)
        print(f"Reference energies 로드 완료: {ref_energies_path}")
    
    # ========================================
    # SECTION 9: CONVERT TO PYMATGEN ENTRIES
    # ========================================
    print(f"\n[진행 중] pymatgen PourbaixEntry 변환 중...")
    
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
            
            entry = PourbaixEntry(ComputedEntry(comp_str, formation_energy), entry_id=data['name'])
            bulk_entries.append(entry)
        except Exception as e:
            print(f"Warning: Failed to create entry for {data['name']} (comp: {comp_str}): {e}")
    
    print(f"  Bulk entries: {len(bulk_entries)}개")
    
    # Solid entries from thermodynamic data
    solid_entries = []
    for el in unique_elements:
        if el in thermo_data and 'solids' in thermo_data[el]:
            for solid_formula, energy in thermo_data[el]['solids'].items():
                try:
                    energy_ev = energy / calmol
                    entry = PourbaixEntry(PDEntry(solid_formula, energy_ev))
                    solid_entries.append(entry)
                except Exception as e:
                    print(f"Warning: Solid '{solid_formula}' 처리 실패: {e}")
    
    print(f"  Solid entries: {len(solid_entries)}개")
    
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
                    print(f"Warning: Ion '{ion_formula}' 처리 실패: {e}")
    
    print(f"  Ion entries: {len(ion_entries)}개")
    
    # Print all entries list
    print("\n" + "="*70)
    print("모든 엔트리 리스트")
    print("="*70)
    
    print(f"\n[Bulk Entries from JSON files] ({len(bulk_entries)}개):")
    for i, entry in enumerate(bulk_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV), entry_id: {entry.entry_id}")
    
    print(f"\n[Solid Entries from thermodynamic_data.json] ({len(solid_entries)}개):")
    for i, entry in enumerate(solid_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV)")
    
    print(f"\n[Ion Entries from thermodynamic_data.json] ({len(ion_entries)}개):")
    for i, entry in enumerate(ion_entries, 1):
        energy = entry.energy if hasattr(entry, 'energy') else entry.entry.energy
        conc = entry.concentration if hasattr(entry, 'concentration') else 'N/A'
        print(f"  {i}. {entry.name} (composition: {entry.composition}, energy: {energy:.6f} eV, conc: {conc})")
    
    # all_entries = solid_entries + ion_entries
    all_entries = solid_entries + bulk_entries + ion_entries
    
    # Remove duplicate entry names, keep only the one with lowest energy
    print(f"\n[진행 중] 중복 entry_name 제거 중...")
    print(f"  원본 엔트리 수: {len(all_entries)}개")
    
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
            print(f"\n  [DEBUG] 중복 발견: entry_name='{entry_name}' ({len(entries_list)}개)")
            for i, (entry, energy, entry_id, comp_str) in enumerate(entries_list):
                status = "✓ 유지" if i == 0 else "✗ 제거"
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
    print(f"\n  제거된 엔트리 수: {removed_count}개")
    print(f"  남은 엔트리 수: {len(all_entries)}개")
    print(f"\n[Total Entries] ({len(all_entries)}개)")
    print("="*70)
    
    # ========================================
    # SECTION 10: PLOT GENERATION
    # ========================================
    print("\n" + "="*70)
    print("다이어그램 생성 시작")
    print("="*70)
    
    # Combined diagram (all entries)
    print("\n[진행 중] Pourbaix 다이어그램 생성 중...")
    plot_pourbaix(all_entries, f'pourbaix{args.suffix}.png')
    
    print("="*70)
    print("모든 다이어그램 생성 완료!")
    print("="*70)

def plot_pourbaix(entries, png_name):
    """Plot Pourbaix diagram using pymatgen PourbaixPlotter"""
    if len(entries) == 0:
        print(f"  Warning: 엔트리가 없어서 다이어그램을 생성할 수 없습니다.")
        return
    
    print(f"\n[진행 중] {png_name} 다이어그램 생성 시작...")
    print(f"  - 엔트리 수: {len(entries)}개")
    
    print(f"  - PourbaixDiagram 계산 중 (convex hull 계산, 시간이 걸릴 수 있습니다)...")
    start_time = time.time()
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    elapsed = time.time() - start_time
    print(f"  ✓ PourbaixDiagram 계산 완료 ({elapsed:.1f}초 소요)")
    print(f"  - 안정 엔트리 수: {len(pourbaix.stable_entries)}개")
    
    print(f"  - Plotter 초기화 중...")
    plotter = PourbaixPlotter(pourbaix)
    
    fig, ax = plt.subplots(figsize=(7, 6))
    plotter.get_pourbaix_plot(limits=[[args.pHmin, args.pHmax], [args.Umin, args.Umax]], 
                              label_domains=True, label_fontsize=10,
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
    
    entry_name_mapping = {}  # Plain format -> entry_id
    entry_latex_mapping = {}  # Plain format -> LaTeX format (for replacement)
    for entry in entries:
        if hasattr(entry, 'entry_id') and entry.entry_id:
            # Convert entry.name from LaTeX to plain format for matching
            plain_name = remove_latex_formatting(entry.name)
            entry_name_mapping[plain_name] = entry.entry_id
            entry_latex_mapping[plain_name] = entry.name  # Store LaTeX format for replacement

    
    # Prepare colors for different entry types
    # Separate entries into categories based on name patterns
    bulk_entries_list = [e for e in stable_entries if '(s)' not in e.name and '(aq)' not in e.name and '[' not in e.name]
    solid_entries_list = [e for e in stable_entries if '(s)' in e.name]
    ion_entries_list = [e for e in stable_entries if '(aq)' in e.name or '[' in e.name]
    
    # Assign colors using colormaps
    n_bulk = len(bulk_entries_list)
    n_solid = len(solid_entries_list)
    n_ion = len(ion_entries_list)
    
    # Use different colormaps for different entry types
    if n_bulk > 0:
        bulk_colors = [plt.cm.Blues(i) for i in np.linspace(0.2, 0.8, n_bulk)]
    else:
        bulk_colors = []
    
    if n_solid > 0:
        solid_colors = [plt.cm.Oranges(i) for i in np.linspace(0.2, 0.8, n_solid)]
    else:
        solid_colors = []
    
    if n_ion > 0:
        ion_colors = [plt.cm.Purples(i) for i in np.linspace(0.2, 0.8, n_ion)]
    else:
        ion_colors = []
    
    # Create color mapping for each entry
    entry_color_map = {}
    for i, entry in enumerate(bulk_entries_list):
        entry_color_map[entry.name] = bulk_colors[i % len(bulk_colors)] if bulk_colors else 'lightblue'
    for i, entry in enumerate(solid_entries_list):
        entry_color_map[entry.name] = solid_colors[i % len(solid_colors)] if solid_colors else 'lightcoral'
    for i, entry in enumerate(ion_entries_list):
        entry_color_map[entry.name] = ion_colors[i % len(ion_colors)] if ion_colors else 'lavender'
    
    # Fill each domain with color
    for entry in stable_entries:
        try:
            vertices = plotter.domain_vertices(entry)
            if len(vertices) > 0:
                x, y = zip(*vertices)
                color = entry_color_map.get(entry.name, 'lightgray')
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
    
    # Add OER line if bulk diagram
    if 'bulk' in png_name:
        pH2 = np.arange(args.pHmin, args.pHmax + 0.01, 0.01)
        plt.plot(pH2, 1.23 - pH2 * const, color='black', lw=2.0)
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    print(f"  - 플롯 저장 중...")
    plt.tight_layout()
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    if args.show_fig:
        plt.show()
    plt.close()
    
    total_time = time.time() - start_time
    print(f"  ✓ {png_name} 다이어그램 생성 완료! (총 {total_time:.1f}초 소요)\n")

if __name__ == "__main__":
    main()
