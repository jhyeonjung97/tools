import os
import re
import socket
import argparse
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
import glob
from ase.io import read
from mendeleev import element
from collections import Counter
import pandas as pd
from pymatgen.core.ion import Ion
import json

# Global constants definition
SUBSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}
SUPERSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}
Diverging_colors = ['RdBu', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--json-dir', type=str, default='.', help='Folder path containing json files (default: current folder)')
parser.add_argument('--csv-dir', type=str, default='.', help='Folder path containing label.csv (default: current folder)')
parser.add_argument('--suffix', type=str, default='', help='Suffix for the output filename (e.g., --suffix "Co")')
parser.add_argument('--hybrid', action='store_true', help='hybrid mode')
parser.add_argument('--gc', action='store_true', help='Apply Grand Canonical DFT using A, B, C columns')
parser.add_argument('--pH', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--conc', type=float, default=10**-6, help='concentration (default: 10^-6)')
parser.add_argument('--gibbs', action='store_true', help='Apply Gibbs free energy correction from G_corr column')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--pHmin', type=float, default=0, help='Minimum pH value (default: 0)')
parser.add_argument('--pHmax', type=float, default=14, help='Maximum pH value (default: 14)')
parser.add_argument('--Umin', type=float, default=-1, help='Minimum potential value (default: -1)')
parser.add_argument('--Umax', type=float, default=3, help='Maximum potential value (default: 3)')
parser.add_argument('--Gmin', type=float, help='Minimum Gibbs free energy value for 1D plot')
parser.add_argument('--Gmax', type=float, help='Maximum Gibbs free energy value for 1D plot')
parser.add_argument('--figx', type=float, default=4, help='Figure width in inches (default: 4)')
parser.add_argument('--figy', type=float, default=3, help='Figure height in inches (default: 3)')
parser.add_argument('--HER', action='store_true', help='Show HER line (0-pH*const)')
parser.add_argument('--OER', action='store_true', help='Show OER line (1.23-pH*const)')
parser.add_argument('--legend-in', action='store_true', help='Show legend inside the plot')
parser.add_argument('--legend-out', action='store_true', help='Show legend outside the plot')
parser.add_argument('--legend-up', action='store_true', help='Show legend above the plot')
parser.add_argument('--cmap', type=str, default='Greys', help='Color map for bulk plot (default: Greys)')
parser.add_argument('--cmin', type=float, default=0.1, help='Minimum value for color range for bulk plot (default: 0.1)')
parser.add_argument('--cmax', type=float, default=0.7, help='Maximum value for color range for bulk plot (default: 0.7)')
parser.add_argument('--cgap', type=float, default=0.0, help='Fraction of colormap to skip in the center for bulk plot (e.g., 0.2 means skip 0.4~0.6)')
parser.add_argument('--cmap-2d', type=str, default='RdBu', help='Color map for 2D plot (default: RdBu)')
parser.add_argument('--cmin-2d', type=float, default=0.0, help='Minimum value for color range for 2D plot (default: 0.0)')
parser.add_argument('--cmax-2d', type=float, default=1.0, help='Maximum value for color range for 2D plot (default: 1.0)')
parser.add_argument('--cgap-2d', type=float, default=0.2, help='Fraction of colormap to skip in the center for 2D plot (e.g., 0.2 means skip 0.4~0.6)')
parser.add_argument('--cmap-1d', type=str, default='Spectral', help='Color map for 1D plot (default: Spectral)')
parser.add_argument('--cmin-1d', type=float, default=0.0, help='Minimum value for color range for 1D plot (default: 0.0)')
parser.add_argument('--cmax-1d', type=float, default=1.0, help='Maximum value for color range for 1D plot (default: 1.0)')
parser.add_argument('--cgap-1d', type=float, default=0.0, help='Fraction of colormap to skip in the center for 1D plot (e.g., 0.2 means skip 0.4~0.6)')
parser.add_argument('--no-bulk', action='store_true', help='Do not show bulk plot')
parser.add_argument('--show-fig', action='store_true', help='Show the plot')
parser.add_argument('--show-thermo', action='store_true', help='Show thermodynamic data details')
parser.add_argument('--show-ref', action='store_true', help='Show reference surface information')
parser.add_argument('--show-element', action='store_true', help='Show elements information')
parser.add_argument('--show-count', action='store_true', help='Show minimum count of each element')
parser.add_argument('--show-label', action='store_true', help='Show file labels')
parser.add_argument('--show-min-coord', action='store_true', help='Show minimum coordinate for each surface')
args = parser.parse_args()

is_hybrid = args.hybrid
is_gibbs = args.gibbs
is_gc = args.gc

png_name = 'pourbaix'
suffix = ''

if hasattr(args, 'hybrid') and is_hybrid:
    png_name += '_hybrid'
else:
    png_name += '_surface'

if hasattr(args, 'gc') and is_gc:
    suffix += '_gc'

if hasattr(args, 'legend_in') and args.legend_in:
    suffix += '_legend_in'
elif hasattr(args, 'legend_out') and args.legend_out:
    suffix += '_legend_out'
elif hasattr(args, 'legend_up') and args.legend_up:
    suffix += '_legend_up'

if hasattr(args, 'suffix') and args.suffix:
    suffix += '_' + args.suffix

def main():
    # Set pH and potential grid
    tick = args.tick if hasattr(args, 'tick') else 0.01
    pHmin, pHmax = args.pHmin, args.pHmax
    pHrange = np.arange(pHmin, pHmax + tick, tick)
    Umin, Umax = args.Umin, args.Umax
    Urange = np.arange(Umin, Umax + 0.06 * 14, tick)
    target_pH = args.pH if hasattr(args, 'pH') else 0
    
    elements = set()
    json_dir = args.json_dir
    json_files = glob.glob(os.path.join(json_dir, "*.json"))
    json_files = sorted(json_files, key=lambda x: os.path.basename(x))

    csv_dir = args.csv_dir
    label_csv_path = os.path.join(csv_dir, 'label.csv')
    
    file_labels = {}
    file_oh_counts = {}
    file_gibbs_corrections = {}
    file_gc_params = {}  # A, B, C
    
    if os.path.exists(label_csv_path):
        # Read CSV file with header
        label_df = pd.read_csv(label_csv_path, header=None, names=['json_name', 'label', '#OH', 'G_corr', 'A', 'B', 'C'])
        for idx, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']
            if '#OH' in row:
                file_oh_counts[json_name] = float(row['#OH'])
            if is_gibbs and 'G_corr' in row:
                file_gibbs_corrections[json_name] = float(row['G_corr'])
            if is_gc and all(col in row for col in ['A', 'B', 'C']):
                file_gc_params[json_name] = {
                    'A': float(row['A']),
                    'B': float(row['B']),
                    'C': float(row['C'])
                }
    else:
        for json_file in json_files:
            atoms = read(json_file)
            formula = atoms.get_chemical_formula()
            file_labels[os.path.basename(json_file)] = formula
    # Sort by atomic number
    for json_file in json_files:
        atoms = read(json_file)
        elements.update(atoms.get_chemical_symbols())
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    if args.show_label:
        print("\nFile labels:")
        for fname, label in file_labels.items():
            print(f"{fname}: {label}")

    # Generate surfs DataFrame
    # List of all elements sorted by atomic number
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    # Dictionary to store minimum count for each element (initial value is None)
    min_counts = {el: None for el in sorted_elements}
    surfs = []
    for json_file in json_files:
        atoms = read(json_file)
        energy = atoms.get_potential_energy()
        symbols = atoms.get_chemical_symbols()
        row = {}
        row['E_DFT'] = energy
        row['e'] = 0
        # Count number of each element  
        symbol_count = Counter(symbols)
        for el in sorted_elements:
            count = symbol_count[el]
            row[el] = float(count)
            # Update minimum value
            if min_counts[el] is None:
                min_counts[el] = count
            else:
                min_counts[el] = min(min_counts[el], count)
        # Store filename or label in name
        row['conc'] = 1
        json_basename = os.path.basename(json_file)
        row['name'] = file_labels.get(json_basename, json_file)
        
        # Add Gibbs correction
        if is_gibbs and json_basename in file_gibbs_corrections:
            row['gibbs_corr'] = file_gibbs_corrections[json_basename]
        else:
            row['gibbs_corr'] = 0.0
            
        # Add GC parameters
        if is_gc and json_basename in file_gc_params:
            row['A'] = file_gc_params[json_basename]['A']
            row['B'] = file_gc_params[json_basename]['B']
            row['E_DFT'] = file_gc_params[json_basename]['C']
        else:
            row['A'] = 0.0
            row['B'] = 0.0    
        surfs.append(row)

    if args.show_count:
        print("\nMinimum count of each element across all files:")
        for el in sorted_elements:
            print(f"{el}: {min_counts[el]}")

    # Subtract minimum count from each row
    for row in surfs:
        for el in sorted_elements:
            row[el] = row[el] - min_counts[el]
    
    # Remove element columns with all values being 0
    zero_cols = [el for el in sorted_elements if all(row[el] == 0 for row in surfs)]
    remaining_elements = [el for el in sorted_elements if el not in zero_cols]
    
    # Delete keys corresponding to zero_cols from surfs
    for row in surfs:
        for col in zero_cols:
            del row[col]
    
    # Define unique_elements as elements excluding O and H
    unique_elements = [el for el in remaining_elements if el not in ['O', 'H']]
    if args.show_element:
        print("\nsorted_elements:", sorted_elements)
        print("remaining_elements:", remaining_elements)
        print("unique_elements:", unique_elements)

    # Set ref_surf as the one with lowest energy among those with all unique_elements = 0
    ref_candidates = []
    for i, row in enumerate(surfs):
        # Check if all elements in unique_elements are 0
        if all(row[elem] == 0.0 for elem in unique_elements if elem in row):
            ref_candidates.append((i, float(row['E_DFT'])))

    if ref_candidates:
        # Select the one with highest energy as ref_surf
        ref_surf_idx, ref_energy = max(ref_candidates, key=lambda x: x[1])
        ref_surf = surfs[ref_surf_idx]
        if args.show_ref:
            print(f"\nReference surface (index {ref_surf_idx}): {ref_surf['name']}, E_DFT: {ref_energy:.2f} eV")
            print(f"Composition: " + ", ".join([f"{elem}: {ref_surf[elem]}" for elem in remaining_elements if elem in ref_surf]))
    else:
        if args.show_ref:
            print("No reference surface found (no surface with all unique_elements = 0)")

    # Apply formation energy correction
    reference_surface_energy = ref_surf['E_DFT']
    for k in range(len(surfs)):
        # Get OH count
        json_basename = None
        for fname, label in file_labels.items():
            if label == surfs[k]['name']:
                json_basename = fname
                break
        
        if json_basename and json_basename in file_oh_counts:
            if not np.isnan(file_oh_counts[json_basename]):
                oh_count = file_oh_counts[json_basename]
            else:
                oh_count = 0
        
        formation_energy_correction = (
            - (surfs[k]['H'] - oh_count) * (gh - dgh)
            - (surfs[k]['O'] - oh_count) * (go - dgo)
            - oh_count * (goh - dgoh)
        )
        
        # Apply Gibbs free energy correction
        gibbs_correction = surfs[k]['gibbs_corr'] if is_gibbs else 0.0
        surfs[k]['E_DFT'] = surfs[k]['E_DFT'] - reference_surface_energy + formation_energy_correction + gibbs_correction

    # Get data corresponding to unique_elements from thermodynamic_data.json
    thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.json')
    thermo_data = {}
    if os.path.exists(thermo_data_path):
        with open(thermo_data_path, 'r') as f:
            thermo_data = json.load(f)
    
    # Read reference_energies.json
    ref_energies_path = os.path.join(os.path.dirname(__file__), 'reference_energies.json')
    ref_energies = {}
    if os.path.exists(ref_energies_path):
        with open(ref_energies_path, 'r') as f:
            ref_energies = json.load(f)
    
    # Create solids, ions, gases, liquids lists
    ions = []
    solids = []
    gases = []
    liquids = []
    
    # Process thermodynamic data for unique_elements
    for el in unique_elements:
        ions_el = []
        solids_el = []
        gases_el = []
        liquids_el = []
        # Check if element exists in reference energies
        if el not in ref_energies:
            print(f"WARNING: Element '{el}' not found in reference_energies.json! Energy corrections may be inaccurate.")
            
        if el in thermo_data:
            if args.show_thermo:
                print(f"\n{el} thermodynamic data:")
                
            # Process ions
            if 'ions' in thermo_data[el] and thermo_data[el]['ions'] != {}:
                if args.show_thermo:
                    print(f"  Ions reduced dict:")
                for ion_formula, energy in thermo_data[el]['ions'].items():
                    try:
                        reduced_dict = Ion.from_formula(ion_formula).to_reduced_dict
                        if args.show_thermo:
                            print(f"    {ion_formula}: {reduced_dict}, energy: {energy}")
                        
                        # Calculate energy correction
                        energy_correction = 0
                        # Add water energy for each O
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # Add energy for elements in reference_energies
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # Create row (using corrected energy)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = args.conc                        
                        row['name'] = format_name(ion_formula) + '(aq)'
                        row['A'] = 0.0
                        row['B'] = 0.0
                        
                        # Normalize by el count
                        el_count = reduced_dict.get(el, 1)  # Set to 1 if el is not present
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # Display fractions as superscript/subscript
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']

                        ions.append(row)
                        ions_el.append(row)
                    except:
                        if args.show_thermo:
                            print(f"    {ion_formula}: parsing failed, energy: {energy}")

            # Process solids
            if 'solids' in thermo_data[el] and thermo_data[el]['solids'] != {}:
                if args.show_thermo:
                    print(f"  Solids reduced dict:")
                for solid_formula, energy in thermo_data[el]['solids'].items():
                    try:
                        reduced_dict = Ion.from_formula(solid_formula).to_reduced_dict
                        if args.show_thermo:
                            print(f"    {solid_formula}: {reduced_dict}, energy: {energy}")
                        
                        # Calculate energy correction
                        energy_correction = 0
                        # Add water energy for each O
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # Add energy for elements in reference_energies
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # Create row (using corrected energy)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = 1
                        row['name'] = format_name(solid_formula) + '(s)'
                        row['A'] = 0.0
                        row['B'] = 0.0

                        # Normalize by el count
                        el_count = reduced_dict.get(el, 1)  # Set to 1 if el is not present
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # Display fractions as superscript/subscript
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']
                        
                        solids.append(row)
                        solids_el.append(row)
                    except:
                        if args.show_thermo:
                            print(f"    {solid_formula}: parsing failed, energy: {energy}")

            # Process gases
            if 'gases' in thermo_data[el] and thermo_data[el]['gases'] != {}:
                if args.show_thermo:
                    print(f"  Gases reduced dict:")
                for gas_formula, energy in thermo_data[el]['gases'].items():
                    try:
                        reduced_dict = Ion.from_formula(gas_formula).to_reduced_dict
                        if args.show_thermo:
                            print(f"    {gas_formula}: {reduced_dict}, energy: {energy}")
                        
                        # Calculate energy correction
                        energy_correction = 0
                        # Add water energy for each O
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # Add energy for elements in reference_energies
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # Create row (using corrected energy)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = args.conc
                        row['name'] = format_name(gas_formula) + '(g)'
                        row['A'] = 0.0
                        row['B'] = 0.0

                        # Normalize by el count
                        el_count = reduced_dict.get(el, 1)  # Set to 1 if el is not present
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # Display fractions as superscript/subscript
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']

                        gases.append(row)
                        gases_el.append(row)
                    except:
                        if args.show_thermo:
                            print(f"    {gas_formula}: parsing failed, energy: {energy}")

            # Process liquids
            if 'liquids' in thermo_data[el] and thermo_data[el]['liquids'] != {}:
                if args.show_thermo:
                    print(f"  Liquids reduced dict:")
                for liquid_formula, energy in thermo_data[el]['liquids'].items():
                    try:
                        reduced_dict = Ion.from_formula(liquid_formula).to_reduced_dict
                        if args.show_thermo:
                            print(f"    {liquid_formula}: {reduced_dict}, energy: {energy}")
                        
                        # Calculate energy correction
                        energy_correction = 0
                        # Add water energy for each O
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # Add energy for elements in reference_energies
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # Create row (using corrected energy)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = 1
                        row['name'] = format_name(liquid_formula) + '(l)'
                        row['A'] = 0.0
                        row['B'] = 0.0

                        # Normalize by el count
                        el_count = reduced_dict.get(el, 1)  # Set to 1 if el is not present
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # Display fractions as superscript/subscript
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']

                        liquids.append(row)
                        liquids_el.append(row)
                    except:
                        if args.show_thermo:
                            print(f"    {liquid_formula}: parsing failed, energy: {energy}")
            
            if is_hybrid and not args.no_bulk:
                bulks = ions_el + solids_el + gases_el + liquids_el
                nbulks = len(bulks)

                fig, ax = plt.subplots(figsize=(args.figx, args.figy))
                ax.axis([pHmin, pHmax, Umin, Umax])
                ax.set_xlabel('pH')
                ax.set_ylabel('E (V vs. SHE)')
                plt.xticks(np.arange(pHmin, pHmax + 1, 2))

                # Calculate lowest bulks
                lowest_bulks = np.full((len(Urange), len(pHrange)), np.nan)
                pHindex = 0
                for pH in pHrange:
                    Uindex = 0
                    for U in Urange:
                        values = []
                        for k in range(nbulks):
                            value = dg(bulks[k], pH, U, bulks[0])
                            values.append(value)
                        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
                        lowest_bulks[Uindex][pHindex] = sorted_values[0]
                        Uindex += 1
                    pHindex += 1

                unique_ids = []
                for i in range(len(Urange)):
                    for j in range(len(pHrange)):
                        bulk_id = int(lowest_bulks[i, j])
                        if bulk_id not in unique_ids:
                            unique_ids.append(bulk_id)

                colormap = getattr(plt.cm, args.cmap, plt.cm.RdBu)
                n_colors = len(unique_ids)
                if args.cgap > 0 and args.cmap in Diverging_colors:
                    gap = args.cgap
                    left_end = 0.5 - gap/2
                    right_start = 0.5 + gap/2
                    n1 = int(np.ceil(n_colors * left_end))
                    n2 = n_colors - n1
                    left = np.linspace(args.cmin, left_end, n1, endpoint=False)
                    right = np.linspace(right_start, args.cmax, n2, endpoint=True)
                    color_values = np.concatenate([left, right])
                else:
                    color_values = np.linspace(args.cmin, args.cmax, n_colors)
                colors = colormap(color_values)
                cmap = mcolors.ListedColormap(colors)
                bounds = np.arange(n_colors + 1) - 0.5
                norm = mcolors.BoundaryNorm(bounds, cmap.N)
                id_map = {val: idx for idx, val in enumerate(unique_ids)}
                mapped_bulks = np.vectorize(id_map.get)(lowest_bulks)

                # Generate legend (bulk, based on pH=0)
                for idx, bulk_id in enumerate(reversed(unique_ids)):
                    label = bulks[int(bulk_id)]['name']
                    plt.plot([], [], color=colors[len(unique_ids)-1-idx], linewidth=5, label=label)

                # pcolormesh
                pH_grid, U = np.meshgrid(pHrange, Urange)
                plt.pcolormesh(pH_grid, U, mapped_bulks, cmap=cmap, norm=norm)

                # Show water stability region
                if args.OER:
                    plt.plot(pHrange, 1.23-pHrange*const, '--', lw=1, color='mediumblue')
                if args.HER:
                    plt.plot(pHrange, 0-pHrange*const, '--', lw=1, color='mediumblue')

                if args.legend_in:
                    plt.legend(fontsize='small', ncol=1, handlelength=3, edgecolor='black', loc='upper right')
                elif args.legend_out:
                    plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., 
                            fontsize='small', ncol=1, handlelength=3, edgecolor='black')
                elif args.legend_up:
                    plt.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', borderaxespad=0., 
                            fontsize='small', ncol=3, handlelength=3, edgecolor='black')

                plt.savefig(f'pourbaix_bulk_{el}{suffix}.png', dpi=300, bbox_inches='tight')
                print(f"Bulk Pourbaix diagram saved as pourbaix_bulk_{el}{suffix}.png")
                if args.show_fig:
                    plt.tight_layout()
                    plt.show()
                plt.close(fig)
        else:
            print(f"WARNING: Element '{el}' not found in thermodynamic data! This element will be ignored in Pourbaix diagram calculations.")
            if args.show_thermo:
                print(f"\n{el}: No thermodynamic data found")
                
    save_surfs = surfs.copy()
    nsurfs, nions, nsolids, ngases, nliquids = len(surfs), len(ions), len(solids), len(gases), len(liquids)
    if is_gc:
        surfs_df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name', 'A', 'B'])
    else:
        surfs_df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])

    print()
    print(format_df_for_display(surfs_df))
    print(f"Surfs: {nsurfs} entries\n")

    if is_hybrid:
        ions_df = pd.DataFrame(ions, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
        solids_df = pd.DataFrame(solids, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
        gases_df = pd.DataFrame(gases, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
        liquids_df = pd.DataFrame(liquids, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
        print(format_df_for_display(ions_df))
        print(f"Ions: {nions} entries\n")
        print(format_df_for_display(solids_df))
        print(f"Solids: {nsolids} entries\n")
        print(format_df_for_display(gases_df))
        print(f"Gases: {ngases} entries\n")
        print(format_df_for_display(liquids_df))
        print(f"Liquids: {nliquids} entries\n")

    new_surfs = []
    
    if is_hybrid:
        # Process each unique_element
        for unique_elem in unique_elements:
            # Find surfs where the unique_element is 0
            for k in range(nsurfs):
                if surfs[k][unique_elem] == 0:
                    # Combine with compounds having 1 unique_element from ions
                    for ion in ions:
                        if ion[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + ion[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * ion[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + ion[key]
                            new_surfs.append(new_surf)
                    
                    # Combine with compounds having 1 unique_element from solids
                    for solid in solids:
                        if solid[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + solid[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * solid[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + solid[key]
                            new_surfs.append(new_surf)
                    
                    # Combine with compounds having 1 unique_element from gases
                    for gas in gases:
                        if gas[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + gas[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * gas[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + gas[key]
                            new_surfs.append(new_surf)
    else:
        # When is_hybrid is False: use only pure forms of each element
        # Find pure form compounds of each element (1 of the element, 0 of others)
        ref_solids = []
        ref_gases = []
        ref_liquids = []
        
        for elem in unique_elements:
            # Find solids with only 1 of this element and 0 of others
            for solid in solids:
                if solid[elem] == 1 and all(solid[other_elem] == 0 for other_elem in remaining_elements if other_elem != elem):
                    if solid not in ref_solids:
                        ref_solids.append(solid)
            
            # Find gases with only 1 of this element and 0 of others
            for gas in gases:
                if gas[elem] == 1 and all(gas[other_elem] == 0 for other_elem in remaining_elements if other_elem != elem):
                    if gas not in ref_gases:
                        ref_gases.append(gas)
            
            # Find liquids with only 1 of this element and 0 of others
            for liquid in liquids:
                if liquid[elem] == 1 and all(liquid[other_elem] == 0 for other_elem in remaining_elements if other_elem != elem):
                    if liquid not in ref_liquids:
                        ref_liquids.append(liquid)
        
        # Process each unique_element
        for unique_elem in unique_elements:
            # Find surfs where the unique_element is 0
            for k in range(nsurfs):
                if surfs[k][unique_elem] == 0:

                    # Combine with compounds having 1 unique_element from ref_solids
                    for solid in ref_solids:
                        if solid[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + solid[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * solid[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + solid[key]
                            new_surfs.append(new_surf)
                    
                    # Combine with compounds having 1 unique_element from ref_gases
                    for gas in ref_gases:
                        if gas[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + gas[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * gas[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + gas[key]
                            new_surfs.append(new_surf)
                    
                    # Combine with compounds having 1 unique_element from ref_liquids
                    for liquid in ref_liquids:
                        if liquid[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + liquid[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * liquid[key]
                                elif key in ['gibbs_corr', 'A', 'B']:
                                    new_surf[key] = surfs[k][key]  # Keep original surface values
                                else:
                                    new_surf[key] = surfs[k][key] + liquid[key]
                            new_surfs.append(new_surf)

    # Add new_surfs to surfs
    surfs.extend(new_surfs)
    # # Keep only those where not all unique_elements are 0
    unique_surfs = [surf for surf in surfs if not all(surf[elem] == 0 for elem in unique_elements)]
    if len(unique_surfs) > 0:
        surfs = unique_surfs
    nsurfs = len(surfs)
    if is_gc:
        df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name', 'A', 'B'])
    else:
        df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    if is_hybrid:
        print(format_df_for_display(df))
        print(f"After adding combinations: {nsurfs} surfs entries\n")

    # Calculate lowest surfaces
    lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

    pHindex = 0
    for pH in pHrange:
        Uindex = 0
        for U in Urange:
            values = []
            for k in range(nsurfs):
                value = dg(surfs[k], pH, U, surfs[ref_surf_idx])
                values.append(value)
            sorted_values = sorted(range(len(values)), key=lambda k: values[k])
            lowest_surfaces[Uindex][pHindex] = sorted_values[0]
            Uindex += 1
        pHindex += 1

    # Find minimum coordinates
    min_coords = {}
    n_rows, n_cols = lowest_surfaces.shape

    for j in range(n_cols):
        for i in range(n_rows):
            sid = int(lowest_surfaces[i, j])
            x = pHrange[j]
            y = Urange[i]
            if sid not in min_coords:
                min_coords[sid] = (x, y)
            else:
                current_x, current_y = min_coords[sid]
                if x < current_x or (x == current_x and y < current_y):
                    min_coords[sid] = (x, y)

    for sid in sorted(min_coords):
        x, y = min_coords[sid]
        name = surfs[int(sid)]['name']
        if args.show_min_coord:
            print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}, name = {name}")

    # Draw Pourbaix diagram
    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    ax.axis([pHmin, pHmax, Umin, Umax])
    ax.set_xlabel('pH')
    ax.set_ylabel('E (V vs. SHE)')
    # ax.tick_params(right=True, direction="in")
    plt.xticks(np.arange(pHmin, pHmax + 1, 2))

    # Generate unique_ids from 2D lowest_surfaces and sort by energy at pH=0
    # Find all unique surface IDs in diagonal order of first appearance in lowest_surfaces
    unique_ids = []
    for i in range(len(Urange)):
        for j in range(len(pHrange)):
            surf_id = int(lowest_surfaces[i, j])
            if surf_id not in unique_ids:
                unique_ids.append(surf_id)

    # Separate unique_ids into two groups: save_surfs vs new combinations
    save_surfs_ids = []
    new_surfs_ids = []
    
    for surf_id in unique_ids:
        # Check if this surface is in save_surfs
        is_in_save_surfs = any(save_surf['name'] == surfs[surf_id]['name'] for save_surf in save_surfs)
        
        if is_in_save_surfs:
            save_surfs_ids.append(surf_id)
        else:
            new_surfs_ids.append(surf_id)

    # Create colormaps for each group
    # Group 1: save_surfs (original surfaces) - use cmap_2d
    colormap_2d = getattr(plt.cm, args.cmap_2d, plt.cm.RdBu)
    n_save = len(save_surfs_ids)
    if n_save > 0:
        if args.cgap_2d > 0 and args.cmap_2d in Diverging_colors:
            gap = args.cgap_2d
            left_end = 0.5 - gap/2
            right_start = 0.5 + gap/2
            n1 = int(np.ceil(n_save * left_end))
            n2 = n_save - n1
            left = np.linspace(args.cmin_2d, left_end, n1, endpoint=False)
            right = np.linspace(right_start, args.cmax_2d, n2, endpoint=True)
            color_values_save = np.concatenate([left, right])
        else:
            color_values_save = np.linspace(args.cmin_2d, args.cmax_2d, n_save)
        colors_save = colormap_2d(color_values_save)
    else:
        colors_save = []

    # Group 2: new combinations - use cmap
    colormap_new = getattr(plt.cm, args.cmap, plt.cm.Greys)
    n_new = len(new_surfs_ids)
    if n_new > 0:
        if args.cgap > 0 and args.cmap in Diverging_colors:
            gap = args.cgap
            left_end = 0.5 - gap/2
            right_start = 0.5 + gap/2
            n1 = int(np.ceil(n_new * left_end))
            n2 = n_new - n1
            left = np.linspace(args.cmin, left_end, n1, endpoint=False)
            right = np.linspace(right_start, args.cmax, n2, endpoint=True)
            color_values_new = np.concatenate([left, right])
        else:
            color_values_new = np.linspace(args.cmin, args.cmax, n_new)
        colors_new = colormap_new(color_values_new)
    else:
        colors_new = []

    # Combine colors and create colormap
    all_colors = []
    id_map = {}
    
    # Add save_surfs colors
    for i, surf_id in enumerate(save_surfs_ids):
        all_colors.append(colors_save[i])
        id_map[surf_id] = len(all_colors) - 1
    
    # Add new_surfs colors
    for i, surf_id in enumerate(new_surfs_ids):
        all_colors.append(colors_new[i])
        id_map[surf_id] = len(all_colors) - 1

    cmap = mcolors.ListedColormap(all_colors)
    bounds = np.arange(len(all_colors) + 1) - 0.5
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    mapped_surfaces = np.vectorize(id_map.get)(lowest_surfaces)

    # Generate legend (2D, based on pH=0)
    # First show save_surfs (reversed)
    for i, surf_id in enumerate(reversed(save_surfs_ids)):
        label = surfs[int(surf_id)]['name']
        plt.plot([], [], color=colors_save[len(save_surfs_ids)-1-i], linewidth=5, label=label)
    
    # Then show new combinations (reversed)
    for i, surf_id in enumerate(reversed(new_surfs_ids)):
        label = surfs[int(surf_id)]['name']
        plt.plot([], [], color=colors_new[len(new_surfs_ids)-1-i], linewidth=5, label=label)

    # pcolormesh
    pH_grid, U = np.meshgrid(pHrange, Urange)
    plt.pcolormesh(pH_grid, U, mapped_surfaces, cmap=cmap, norm=norm)

    # Show water stability region
    if args.OER:
        plt.plot(pHrange, 1.23-pHrange*const, '--', lw=1, color='mediumblue')
    if args.HER:
        plt.plot(pHrange, 0-pHrange*const, '--', lw=1, color='mediumblue')

    if args.legend_in:
        plt.legend(fontsize='small', ncol=1, handlelength=3, edgecolor='black', loc='upper right')
    elif args.legend_out:
        plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., 
                fontsize='small', ncol=1, handlelength=3, edgecolor='black')
    elif args.legend_up:
        plt.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', borderaxespad=0., 
                fontsize='small', ncol=3, handlelength=3, edgecolor='black')

    plt.savefig(f'{png_name}{suffix}.png', dpi=300, bbox_inches='tight')
    print(f"Pourbaix diagram saved as {png_name}{suffix}.png")
    if args.show_fig:
        plt.tight_layout()
        plt.show()

    # === 1D Pourbaix diagram at specific pH ===
    # Calculate lowest and second lowest surfaces at target_pH for 1D plot
    lowest_surfaces_pH = np.zeros(len(Urange))
    second_lowest_surfaces_pH = np.zeros(len(Urange))
    for Uindex, U in enumerate(Urange):
        values = []
        for k in range(nsurfs):
            value = dg(surfs[k], pH=target_pH, U=U, ref_surf=surfs[ref_surf_idx])
            values.append(value)
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces_pH[Uindex] = sorted_values[0]
        if len(sorted_values) > 1:
            second_lowest_surfaces_pH[Uindex] = sorted_values[1]
    
    # Calculate energy at specific pH
    all_energies = []
    unique_second_ids_pH = np.unique(second_lowest_surfaces_pH.astype(int))
    unique_ids_set = set(unique_ids)
    unique_second_ids_set = set(unique_second_ids_pH)

    # 1D plot: Sort unique_ids_set, unique_second_ids_set by 'e' value first, then by energy at U=0
    energies_at_U0 = [(k,surfs[k]['e'] - surfs[k]['H'] + 2*surfs[k]['O'], dg(surfs[k], pH=target_pH, U=0, ref_surf=surfs[ref_surf_idx])) for k in unique_ids_set | unique_second_ids_set]
    sorted_unique_ids = [k for k, _, _ in sorted(energies_at_U0, key=lambda x: (x[1], x[2]))]  # Sort by 'e' first, then by energy

    fig2, ax2 = plt.subplots(figsize=(args.figx, args.figy))
    ax2.axis([Umin, Umax, None, None])
    ax2.set_xlabel('Potential (V vs. SHE)')
    ax2.set_ylabel('Relative Energy (ΔG, eV)')
    ax2.tick_params()

    Urange = np.arange(Umin, Umax, tick)
    
    # 1D plot colormap/color range
    colormap_1d = getattr(plt.cm, args.cmap_1d, plt.cm.RdBu)
    n_colors_1d = len(sorted_unique_ids)
    if args.cgap_1d > 0 and args.cmap_1d in Diverging_colors:
        gap = args.cgap_1d
        left_end = 0.5 - gap/2
        right_start = 0.5 + gap/2
        n1 = int(np.ceil(n_colors_1d * left_end))
        n2 = n_colors_1d - n1
        left = np.linspace(args.cmin_1d, left_end, n1, endpoint=False)
        right = np.linspace(right_start, args.cmax_1d, n2, endpoint=True)
        color_values_1d = np.concatenate([left, right])
    else:
        color_values_1d = np.linspace(args.cmin_1d, args.cmax_1d, n_colors_1d)
    colors_1d = colormap_1d(color_values_1d)
    
    for k in range(nsurfs):
        energies = np.zeros(len(Urange))
        for i, U in enumerate(Urange):
            value = dg(surfs[k], pH=target_pH, U=U, ref_surf=surfs[ref_surf_idx])
            energies[i] = value
        all_energies.extend(energies)
        if not k in sorted_unique_ids:
            ax2.plot(Urange, energies, color='lightgray', alpha=0.3, lw=0.5)
    
    for j, k in enumerate(sorted_unique_ids):
        energies = np.zeros(len(Urange))
        for i, U in enumerate(Urange):
            value = dg(surfs[k], pH=target_pH, U=U, ref_surf=surfs[ref_surf_idx])
            energies[i] = value
        if k in unique_ids_set:
            ax2.plot(Urange, energies, label=surfs[k]['name'], lw=1, color=colors_1d[j])
        elif k in unique_second_ids_set:
            ax2.plot(Urange, energies, '--', label=surfs[k]['name'], lw=1, color=colors_1d[j])
        else:
            ax2.plot(Urange, energies, color='lightgray', alpha=0.3, lw=0.5)

    # Adjust y-axis range
    if args.Gmin:
        Gmin = args.Gmin
    else:
        Gmin = min(all_energies)
    if args.Gmax:
        Gmax = args.Gmax
    else:
        Gmax = max(all_energies)

    ax2.set_ylim(Gmin, Gmax)
    ax2.set_xlim(Umin, Umax)
    # Generate legend (1D)
    if args.legend_in:
        ax2.legend(fontsize='small', ncol=1, handlelength=3, edgecolor='black', loc='upper right')
    elif args.legend_out:
        ax2.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., 
                fontsize='small', ncol=1, handlelength=3, edgecolor='black')
    elif args.legend_up:
        ax2.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', borderaxespad=0., 
                fontsize='small', ncol=3, handlelength=3, edgecolor='black')
    plt.savefig(f'{png_name}_pH{target_pH}{suffix}.png', dpi=300, bbox_inches='tight')
    print(f"Pourbaix diagram saved as {png_name}_pH{target_pH}{suffix}.png")
    if args.show_fig:
        plt.tight_layout()
        plt.show()

# Define dg function (Gibbs free energy calculation)
def dg(surf, pH, U, ref_surf):
    surface_term = (surf['A']*(U**2) + surf['B']*U + surf['E_DFT']) - (ref_surf['A']*(U**2) + ref_surf['B']*U + ref_surf['E_DFT'])
    U_coeff = surf['H'] - 2*surf['O'] - surf['e']
    pH_coeff = surf['H'] - 2*surf['O']
    dg_value = surface_term + U_coeff*U + const*pH_coeff*pH + const*log10(surf['conc'])
    return dg_value

# Convert ion names by counting + or - to superscript, and numbers to subscript
def format_name(formula):
    # Count + or -
    plus_count = formula.count('+')
    minus_count = formula.count('-')
    
    # Base formula with + or - removed
    base_formula = formula.replace('+', '').replace('-', '')
    
    # Convert numbers to subscript
    formatted_formula = ''
    for char in base_formula:
        if char.isdigit():
            formatted_formula += SUBSCRIPT_NUMS[char]
        else:
            formatted_formula += char
    
    # Add charge display
    if plus_count > 0:
        if plus_count == 1:
            return formatted_formula + '⁺'
        else:
            # Convert charge number to superscript
            superscript_num = SUPERSCRIPT_NUMS.get(str(plus_count), str(plus_count))
            return formatted_formula + superscript_num + '⁺'
    elif minus_count > 0:
        if minus_count == 1:
            return formatted_formula + '⁻'
        else:
            # Convert charge number to superscript
            superscript_num = SUPERSCRIPT_NUMS.get(str(minus_count), str(minus_count))
            return formatted_formula + superscript_num + '⁻'
    else:
        return formatted_formula

# Function to format only conc column with scientific notation
def format_df_for_display(df):
    df_copy = df.copy()
    
    # Format numeric columns to 2 decimal places
    numeric_cols = df_copy.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        if col == 'conc':
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.0e}")
        else:
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.2f}")
    
    return df_copy

if __name__ == "__main__":
    # units
    kjmol = 96.485
    calmol = 23.061

    # constants
    kb = 8.617e-5 # eV/K
    T = 298.15 # K
    const = kb * T * np.log(10) # 0.0592 eV
    water = 56.690/calmol

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

    gh = gh2 / 2
    go = gh2o - gh2
    goh = gh2o - gh2 / 2
    gooh = 2 * gh2o - 1.5 * gh2

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

    main()