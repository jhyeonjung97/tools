#!/usr/bin/env python3
"""
HybridPB: Hybrid Pourbaix Diagram Generation Tool

This script generates Pourbaix diagrams (potential-pH diagrams) for electrochemical systems,
with support for both surface and bulk phases. The tool can handle hybrid calculations
incorporating Grand Canonical DFT corrections and thermodynamic data.

Key Features:
- Surface and bulk Pourbaix diagram generation
- Grand Canonical DFT corrections (--gc flag)
- Thermodynamic data integration
- Customizable visualization options
- Support for aqueous species and gas phases
- Export capabilities for publication-ready figures

Author: Hyeonjung Jung
Last modified: 2025

Copyright (C) 2025 Hyeonjung Jung

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard library imports
import os
import re
import argparse
import json
import glob
from math import log10
from collections import Counter

# Scientific computing imports
import numpy as np
import pandas as pd

# Visualization imports
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Materials science and chemistry imports
from ase.io import read
from mendeleev import element
from pymatgen.core.ion import Ion

def load_jsonc(file_path):
    """
    Load JSON file with comments (JSONC format).
    
    This function reads a JSON file that may contain single-line comments
    starting with '//', which are not supported by standard JSON parsers.
    
    Args:
        file_path (str): Path to the JSONC file
        
    Returns:
        dict: Parsed JSON data with comments removed
        
    Raises:
        json.JSONDecodeError: If the file is not valid JSON after comment removal
        FileNotFoundError: If the specified file doesn't exist
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Remove single-line comments (// ...)
    content = re.sub(r'//.*', '', content)
    
    # Parse the cleaned JSON
    return json.loads(content)

# Global constants for chemical formula formatting
# Unicode subscript characters for chemical formulas
SUBSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}

# Unicode superscript characters for ion charges
SUPERSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}

# Available diverging color schemes for visualization
Diverging_colors = ['RdBu', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

# Command line argument parser configuration
# This extensive argument parser allows users to customize every aspect of the Pourbaix diagram generation
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram for electrochemical systems')

# Input/Output file paths
parser.add_argument('--json-dir', type=str, default='.', 
                    help='Folder path containing JSON structure files (default: current folder)')
parser.add_argument('--csv-dir', type=str, default='.', 
                    help='Folder path containing label.csv with species information (default: current folder)')
parser.add_argument('--suffix', type=str, default='', 
                    help='Suffix for the output filename (e.g., --suffix "Co" for Co_pourbaix.png)')

# Calculation modes
parser.add_argument('--hybrid', action='store_true', 
                    help='Enable hybrid mode for surface-bulk hybrid calculations')
parser.add_argument('--gc', action='store_true', 
                    help='Apply Grand Canonical DFT corrections using A, B, C columns from label.csv')

# Thermodynamic conditions
parser.add_argument('--pH', type=int, default=0, 
                    help='pH value for the plot (default: 0)')
parser.add_argument('--concentration', type=float, default=1e-6, 
                    help='Ion concentration in M (default: 10^-6 M)')
parser.add_argument('--pressure', type=float, default=1e-6, 
                    help='Gas pressure in atm (default: 10^-6 atm)')

# Plot range and resolution
parser.add_argument('--tick', type=float, default=0.01, 
                    help='Grid resolution for pH and potential ranges (default: 0.01)')
parser.add_argument('--pHmin', type=float, default=0, 
                    help='Minimum pH value for the diagram (default: 0)')
parser.add_argument('--pHmax', type=float, default=14, 
                    help='Maximum pH value for the diagram (default: 14)')
parser.add_argument('--Umin', type=float, default=-1, 
                    help='Minimum potential value in V vs SHE (default: -1)')
parser.add_argument('--Umax', type=float, default=3, 
                    help='Maximum potential value in V vs SHE (default: 3)')

# 1D plot specific ranges
parser.add_argument('--Gmin', type=float, 
                    help='Minimum Gibbs free energy value for 1D stability plots')
parser.add_argument('--Gmax', type=float, 
                    help='Maximum Gibbs free energy value for 1D stability plots')

# Figure dimensions and layout
parser.add_argument('--figx', type=float, default=4, 
                    help='Figure width in inches (default: 4)')
parser.add_argument('--figy', type=float, default=4, 
                    help='Figure height in inches (default: 4)')

# Electrochemical reference lines
parser.add_argument('--HER', action='store_true', 
                    help='Show Hydrogen Evolution Reaction line (H+/H2 equilibrium)')
parser.add_argument('--OER', action='store_true', 
                    help='Show Oxygen Evolution Reaction line (O2/H2O equilibrium)')

# Legend positioning options
parser.add_argument('--legend-in', action='store_true', 
                    help='Show legend inside the plot area')
parser.add_argument('--legend-out', action='store_true', 
                    help='Show legend outside the plot area')
parser.add_argument('--legend-up', action='store_true', 
                    help='Show legend above the plot')

# Color scheme customization for bulk plots
parser.add_argument('--cmap', type=str, default='Greys', 
                    help='Matplotlib colormap for bulk stability regions (default: Greys)')
parser.add_argument('--cmin', type=float, default=0.1, 
                    help='Minimum value for color mapping in bulk plots (default: 0.1)')
parser.add_argument('--cmax', type=float, default=0.7, 
                    help='Maximum value for color mapping in bulk plots (default: 0.7)')
parser.add_argument('--cgap', type=float, default=0.0, 
                    help='Fraction of colormap to skip in the center (e.g., 0.2 skips middle 20 percent)')

# Color scheme customization for 2D plots
parser.add_argument('--cmap-2d', type=str, default='RdBu', 
                    help='Matplotlib colormap for 2D phase diagrams (default: RdBu)')
parser.add_argument('--cmin-2d', type=float, default=0.0, 
                    help='Minimum value for color mapping in 2D plots (default: 0.0)')
parser.add_argument('--cmax-2d', type=float, default=1.0, 
                    help='Maximum value for color mapping in 2D plots (default: 1.0)')
parser.add_argument('--cgap-2d', type=float, default=0.2, 
                    help='Color gap for 2D plots (default: 0.2)')

# Color scheme customization for 1D plots
parser.add_argument('--cmap-1d', type=str, default='Spectral', 
                    help='Matplotlib colormap for 1D stability plots (default: Spectral)')
parser.add_argument('--cmin-1d', type=float, default=0.0, 
                    help='Minimum value for color mapping in 1D plots (default: 0.0)')
parser.add_argument('--cmax-1d', type=float, default=1.0, 
                    help='Maximum value for color mapping in 1D plots (default: 1.0)')
parser.add_argument('--cgap-1d', type=float, default=0.0, 
                    help='Color gap for 1D plots (default: 0.0)')

# Display options
parser.add_argument('--no-bulk', action='store_true', 
                    help='Skip bulk phase calculations and display only surface phases')
parser.add_argument('--show-fig', action='store_true', 
                    help='Display the generated plot in a window')

# Debug and information display options
parser.add_argument('--show-thermo', action='store_true', 
                    help='Print detailed thermodynamic data information')
parser.add_argument('--show-ref', action='store_true', 
                    help='Print reference surface selection information')
parser.add_argument('--show-element', action='store_true', 
                    help='Print element composition analysis')
parser.add_argument('--show-count', action='store_true', 
                    help='Print minimum count of each element across all structures')
parser.add_argument('--show-label', action='store_true', 
                    help='Print file labels and their assignments')
parser.add_argument('--show-min-coord', action='store_true', 
                    help='Print minimum coordination analysis for each surface')
parser.add_argument('--show-transitions', action='store_true', 
                    help='Print most stable surface transition points in 1D plots')

# Structure visualization
parser.add_argument('--png', action='store_true', 
                    help='Generate PNG images of all JSON structure files')
parser.add_argument('--png-rotation', type=str, default='-90x, -90y, 0z', 
                    help='Rotation angles for PNG structure images (default: "-90x, -90y, 0z")')

# Custom file paths (override defaults)
parser.add_argument('--label-csv', type=str, 
                    help='Custom path to label.csv file (default: ./label.csv)')
parser.add_argument('--thermo-data', type=str, 
                    help='Custom path to thermodynamic data file (default: ./thermodynamic_data.json)')
parser.add_argument('--ref-energies', type=str, 
                    help='Custom path to reference energies file (default: ./reference_energies.json)')

# Parse command line arguments
args = parser.parse_args()

# Extract key calculation modes from arguments
is_hybrid = args.hybrid      # Hybrid surface-bulk calculations
is_gc = args.gc              # Apply Grand Canonical DFT corrections

# Build output filename based on calculation modes and options
png_name = 'pourbaix'
suffix = ''

# Add mode indicators to filename
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
    """
    Main function that orchestrates the entire Pourbaix diagram generation process.
    
    This function handles:
    1. Grid setup for pH and potential ranges
    2. File discovery and structure reading
    3. Label processing from CSV files
    4. Energy calculations and corrections
    5. Thermodynamic data integration
    6. Bulk and surface phase analysis
    7. Plot generation and export
    """
    
    # ========================================
    # SECTION 1: GRID SETUP AND INITIALIZATION
    # ========================================
    
    # Set up pH and potential grid for diagram calculation
    # The grid resolution is controlled by the --tick argument
    tick = args.tick if hasattr(args, 'tick') else 0.01
    pHmin, pHmax = args.pHmin, args.pHmax
    pHrange = np.arange(pHmin, pHmax + tick, tick)  # pH range array
    
    # Potential range setup (V vs SHE - Standard Hydrogen Electrode)
    Umin, Umax = args.Umin, args.Umax
    # Extra range added (0.06 * 14) to accommodate thermodynamic boundaries
    Urange = np.arange(Umin, Umax + 0.06 * 14, tick)  # Potential range array
    
    # Target pH for specific calculations
    target_pH = args.pH if hasattr(args, 'pH') else 0
    
    # ========================================
    # SECTION 2: FILE DISCOVERY AND STRUCTURE READING
    # ========================================
    
    # Initialize element tracking and file discovery
    elements = set()  # Track all unique elements found in structures
    json_dir = args.json_dir
    
    # Find and sort all JSON structure files
    json_files = glob.glob(os.path.join(json_dir, "*.json"))
    json_files = sorted(json_files, key=lambda x: os.path.basename(x))

    # ========================================
    # SECTION 3: LABEL PROCESSING FROM CSV FILES
    # ========================================
    
    # Determine path to label.csv file containing species information
    csv_dir = args.csv_dir
    if hasattr(args, 'label_csv') and args.label_csv:
        label_csv_path = args.label_csv  # Use custom path if provided
    else:
        label_csv_path = os.path.join(csv_dir, 'label.csv')  # Default path
    
    # Initialize dictionaries to store species information
    file_labels = {}              # Maps filenames to display labels
    file_oh_counts = {}           # Maps filenames to OH group counts
    file_gibbs_corrections = {}   # Maps filenames to Gibbs free energy corrections
    file_gc_params = {}           # Maps filenames to Grand Canonical DFT parameters (A, B, C)
    
    # Process label.csv file if it exists
    if os.path.exists(label_csv_path):
        # Read CSV file with predefined column structure:
        # json_name, label, #OH, G_corr, A, B, C
        label_df = pd.read_csv(label_csv_path, header=None, 
                              names=['json_name', 'label', '#OH', 'G_corr', 'A', 'B', 'C'])
        
        # Process each row in the CSV file
        for idx, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']  # Store display label
            
            # Store OH group count if available
            if '#OH' in row and not pd.isna(row['#OH']):
                file_oh_counts[json_name] = float(row['#OH'])
            
            # Store Gibbs free energy correction if available
            if 'G_corr' in row and not pd.isna(row['G_corr']):
                file_gibbs_corrections[json_name] = float(row['G_corr'])
            
            # Store Grand Canonical DFT parameters if GC mode is enabled
            if is_gc and all(col in row for col in ['A', 'B', 'C']):
                if all(not pd.isna(row[col]) for col in ['A', 'B', 'C']):
                    file_gc_params[json_name] = {
                        'A': float(row['A']),  # Chemical potential coefficient
                        'B': float(row['B']),  # pH-dependent term
                        'C': float(row['C']),  # Energy correction
                    }
    else:
        # Fallback: Generate labels from chemical formulas if no CSV file exists
        print(f"Warning: Label file {label_csv_path} not found. Using chemical formulas as labels.")
        for json_file in json_files:
            atoms = read(json_file)
            formula = atoms.get_chemical_formula()
            file_labels[os.path.basename(json_file)] = formula
    # ========================================
    # SECTION 4: ELEMENT ANALYSIS AND STRUCTURE VISUALIZATION
    # ========================================
    
    # Analyze all elements present in the structures
    for json_file in json_files:
        atoms = read(json_file)
        elements.update(atoms.get_chemical_symbols())  # Collect unique elements
    
    # Sort elements by atomic number for consistent ordering
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    
    # Display file labels if requested (debug option)
    if args.show_label:
        print("\nFile labels:")
        for fname, label in file_labels.items():
            print(f"{fname}: {label}")

    # Generate PNG images of structures if requested
    if args.png:
        print("\nSaving PNG images of structures...")
        for json_file in json_files:
            atoms = read(json_file)
            json_basename = os.path.basename(json_file)
            png_filename = json_basename.replace('.json', '.png')
            try:
                from ase.io import write
                # Generate PNG with specified rotation and no unit cell display
                write(png_filename, atoms, rotation=args.png_rotation, show_unit_cell=0)
                print(f"  Generated: {png_filename}")
            except Exception as e:
                print(f"  Failed to save {png_filename}: {e}")
        print("PNG image generation completed.\n")

    # ========================================
    # SECTION 5: SURFACE DATA STRUCTURE CONSTRUCTION
    # ========================================
    
    # Build comprehensive surface data structure
    # This section creates the main data structure that will be used for all calculations
    
    # Re-sort elements by atomic number for consistency
    # Ensure H and O are always included for electrochemical calculations
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    if 'H' not in sorted_elements:
        sorted_elements.append('H')
    if 'O' not in sorted_elements:
        sorted_elements.append('O')
    # Re-sort to maintain atomic number ordering
    sorted_elements = sorted(sorted_elements, key=lambda x: element(x).atomic_number)
    
    # Track minimum count of each element across all structures
    # This is used later for reference energy calculations
    min_counts = {el: None for el in sorted_elements}
    
    # Main surface data list - each entry represents one structure
    surfs = []
    
    print(f"Processing {len(json_files)} structure files...")
    
    for json_file in json_files:
        # Read structure and extract basic information
        atoms = read(json_file)
        energy = atoms.get_potential_energy()  # DFT total energy
        symbols = atoms.get_chemical_symbols()  # List of atomic symbols
        
        # Initialize data row for this structure
        row = {}
        row['E_DFT'] = energy  # DFT total energy (eV)
        row['e'] = 0           # Number of electrons (will be calculated later)
        
        # Count atoms of each element in this structure
        symbol_count = Counter(symbols)
        for el in sorted_elements:
            count = symbol_count.get(el, 0)  # Get count, default to 0 if element not present
            row[el] = float(count)
            
            # Track minimum count across all structures for each element
            if min_counts[el] is None:
                min_counts[el] = count
            else:
                min_counts[el] = min(min_counts[el], count)
        
        # Store metadata
        row['conc'] = 1  # Concentration (default: 1, may be modified later)
        json_basename = os.path.basename(json_file)
        row['name'] = file_labels.get(json_basename, json_file)  # Display name
            
        # Apply Grand Canonical DFT parameters if available
        if is_gc and json_basename in file_gc_params:
            # GC-DFT modifies the energy based on chemical potentials
            gc_params = file_gc_params[json_basename]
            row['A'] = gc_params['A']  # Chemical potential coefficient
            row['B'] = gc_params['B']  # pH-dependent correction
            row['E_DFT'] = gc_params['C']  # Use corrected energy instead of DFT energy
        else:
            row['A'] = 0.0
            row['B'] = 0.0 

        # Apply Gibbs free energy correction if available
        if json_basename in file_gibbs_corrections:
            row['E_DFT'] += file_gibbs_corrections[json_basename]

        surfs.append(row)

    # ========================================
    # SECTION 6: ELEMENT COUNT NORMALIZATION
    # ========================================
    
    # Display minimum element counts if requested (debug option)
    if args.show_count:
        print("\nMinimum count of each element across all files:")
        for el in sorted_elements:
            print(f"{el}: {min_counts[el]}")

    # Normalize element counts by subtracting minimum counts
    # This step converts absolute counts to relative counts for formation energy calculations
    if args.show_count:
        print("Normalizing element counts relative to reference...")
    for row in surfs:
        for el in sorted_elements:
            row[el] = row[el] - min_counts[el]
    
    # Remove element columns that are zero for all structures after normalization
    # These elements don't contribute to the phase equilibria
    # Note: H and O are always kept for electrochemical calculations even if count is 0
    zero_cols = [el for el in sorted_elements if el not in ['H', 'O'] and all(row[el] == 0 for row in surfs)]
    remaining_elements = [el for el in sorted_elements if el not in zero_cols]
    
    if zero_cols and args.show_element:
        print(f"Removing elements with zero variation: {zero_cols}")
    
    # Clean up the data structure by removing zero columns
    for row in surfs:
        for col in zero_cols:
            del row[col]
    
    # Define unique_elements as non-H, non-O elements for phase diagram construction
    # H and O are typically treated as reservoir species in electrochemical systems
    unique_elements = [el for el in remaining_elements if el not in ['O', 'H']]
    
    # Display element analysis if requested (debug option)
    if args.show_element:
        print(f"\nElement analysis:")
        print(f"  All elements found: {sorted_elements}")
        print(f"  Elements with variation: {remaining_elements}")  
        print(f"  Phase-determining elements (non-H, non-O): {unique_elements}")

    # ========================================
    # SECTION 7: REFERENCE SURFACE SELECTION
    # ========================================
    
    # Select reference surface for formation energy calculations
    # The reference surface should have all unique_elements = 0 (i.e., contains only H and O)
    # Among such surfaces, we choose the one with the highest energy for numerical stability
    
    ref_candidates = []
    if args.show_ref:
        print("Searching for reference surface candidates...")
    
    for i, row in enumerate(surfs):
        # Check if this surface contains only H and O (all unique elements = 0)
        if all(row.get(elem, 0) == 0.0 for elem in unique_elements):
            ref_candidates.append((i, float(row['E_DFT']), row['name']))
            if args.show_ref:
                print(f"  Candidate {i}: {row['name']} (E_DFT: {row['E_DFT']:.3f} eV)")

    if ref_candidates:
        # Select the candidate with highest energy as reference surface
        # Higher energy reference makes all formation energies more negative, improving numerical stability
        ref_surf_idx, ref_energy, ref_name = max(ref_candidates, key=lambda x: x[1])
        ref_surf = surfs[ref_surf_idx]
        
        if args.show_ref:
            print(f"\nSelected reference surface: {ref_name} (index {ref_surf_idx})")
            print(f"  Reference energy: {ref_energy:.3f} eV")
        
        if args.show_ref:
            print(f"  Detailed composition:")
            for elem in remaining_elements:
                if elem in ref_surf:
                    print(f"    {elem}: {ref_surf[elem]}")
    else:
        # No suitable reference surface found - this is problematic
        ref_surf = None
        print("WARNING: No reference surface found!")
        print("  A reference surface should contain only H and O atoms.")
        print("  This may indicate an issue with the input structures.")
        
        if args.show_ref:
            print("\nAll surface compositions:")
            for i, row in enumerate(surfs):
                composition = ", ".join([f"{elem}: {row[elem]}" for elem in unique_elements if elem in row and row[elem] != 0])
                print(f"  {i}: {row['name']} - {composition if composition else 'Only H/O'}")

    # ========================================
    # SECTION 8: FORMATION ENERGY CORRECTIONS
    # ========================================
    
    # Apply formation energy corrections to convert DFT total energies to formation energies
    # This is crucial for comparing structures with different compositions
    
    if ref_surf is None:
        print("ERROR: Cannot proceed without a reference surface!")
        return
    
    reference_surface_energy = ref_surf['E_DFT']
    if args.show_ref:
        print(f"\nApplying formation energy corrections...")
        print(f"Using reference energy: {reference_surface_energy:.3f} eV")
    
    for k in range(len(surfs)):
        if not file_gibbs_corrections or all(pd.isna(value) for value in file_gibbs_corrections.values()):  # If no Gibbs corrections are provided or all values are NaN
            # Determine OH group count for this surface
            oh_count = 0  # Default: no OH groups
            json_basename = None
            
            # Find the corresponding filename for this surface
            for fname, label in file_labels.items():
                if label == surfs[k]['name']:
                    json_basename = fname
                    break
            
            # Use OH count from CSV file if available
            if json_basename and json_basename in file_oh_counts:
                if not np.isnan(file_oh_counts[json_basename]):
                    oh_count = file_oh_counts[json_basename]
            
            # Calculate formation energy correction using thermodynamic reference states
            # The correction accounts for the chemical potentials of H, O, and OH
            formation_energy_correction = (
                - (surfs[k]['H'] - oh_count) * (gh - dgh)  # H atoms (excluding those in OH)
                - (surfs[k]['O'] - oh_count) * (go - dgo)  # O atoms (excluding those in OH) 
                - oh_count * (goh - dgoh)                  # OH groups
            )
        else:  # If Gibbs corrections are provided, use simpler correction
            formation_energy_correction = (
                - surfs[k]['H'] * gh  # H atoms
                - surfs[k]['O'] * go  # O atoms
            )
        
        # Update the energy to formation energy
        original_energy = surfs[k]['E_DFT']
        surfs[k]['E_DFT'] = (original_energy - reference_surface_energy + 
                           formation_energy_correction)
        
        if args.show_ref and k < 5:  # Show details for first few surfaces
            print(f"  {surfs[k]['name']}: {original_energy:.3f} → {surfs[k]['E_DFT']:.3f} eV")

    # ========================================
    # SECTION 9: THERMODYNAMIC DATA LOADING
    # ========================================
    
    # Load thermodynamic data for bulk species (ions, gases, solids)
    # This data is essential for constructing the complete Pourbaix diagram
    
    # Determine path to thermodynamic data file
    if hasattr(args, 'thermo_data') and args.thermo_data:
        thermo_data_path = args.thermo_data
    else:
        # Default: look in same directory as script
        thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.json')
        # Fallback to .jsonc extension
        if not os.path.exists(thermo_data_path):
            thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.jsonc')
    
    if args.show_thermo:
        print(f"\nLoading thermodynamic data from: {thermo_data_path}")
    thermo_data = {}
    if os.path.exists(thermo_data_path):
        thermo_data = load_jsonc(thermo_data_path)
    
    # Read reference_energies.json
    ref_energies_path = args.ref_energies if hasattr(args, 'ref_energies') and args.ref_energies else os.path.join(os.path.dirname(__file__), 'reference_energies.json')
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
                        row['conc'] = args.concentration                        
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
                        row['conc'] = args.pressure
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
                ax.set_ylabel('E (V vs. RHE)')
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

                colormap = getattr(plt.cm, args.cmap, plt.cm.Greys)
                n_colors = len(unique_ids)
                if args.cgap > 0 and args.cmap in Diverging_colors:
                    gap = args.cgap
                    left_end = 0.5 - gap/2
                    right_start = 0.5 + gap/2
                    
                    # Calculate how many colors go to left and right
                    n_left = int(np.ceil(n_colors / 2))
                    n_right = n_colors - n_left
                    
                    # Create color values avoiding the gap
                    if n_left > 0:
                        left = np.linspace(args.cmin, left_end, n_left, endpoint=True)
                    else:
                        left = []
                    
                    if n_right > 0:
                        right = np.linspace(right_start, args.cmax, n_right, endpoint=True)
                    else:
                        right = []
                    
                    color_values = np.concatenate([left, right]) if len(left) > 0 and len(right) > 0 else (left if len(left) > 0 else right)
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
    n_original_surfs = len(surfs)  # Store the number of original surfaces
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
                                elif key in ['A', 'B']:
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
                                elif key in ['A', 'B']:
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
                                elif key in ['A', 'B']:
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
                                elif key in ['A', 'B']:
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
                                elif key in ['A', 'B']:
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
                                elif key in ['A', 'B']:
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
    if is_hybrid and args.show_thermo:
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
    ax.set_ylabel('E (V vs. RHE)')
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
        surf_name = surfs[surf_id]['name']
        is_in_save_surfs = any(save_surf['name'] == surf_name for save_surf in save_surfs)
        
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
            
            # Calculate how many colors go to left and right
            n_left = int(np.ceil(n_save / 2))
            n_right = n_save - n_left
            
            # Create color values avoiding the gap
            if n_left > 0:
                left = np.linspace(args.cmin_2d, left_end, n_left, endpoint=True)
            else:
                left = []
            
            if n_right > 0:
                right = np.linspace(right_start, args.cmax_2d, n_right, endpoint=True)
            else:
                right = []
            
            color_values_save = np.concatenate([left, right]) if len(left) > 0 and len(right) > 0 else (left if len(left) > 0 else right)
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
            
            # Calculate how many colors go to left and right
            n_left = int(np.ceil(n_new / 2))
            n_right = n_new - n_left
            
            # Create color values avoiding the gap
            if n_left > 0:
                left = np.linspace(args.cmin, left_end, n_left, endpoint=True)
            else:
                left = []
            
            if n_right > 0:
                right = np.linspace(right_start, args.cmax, n_right, endpoint=True)
            else:
                right = []
            
            color_values_new = np.concatenate([left, right]) if len(left) > 0 and len(right) > 0 else (left if len(left) > 0 else right)
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

    # Find transition points where most stable surface changes (if requested)
    if args.show_transitions:
        print(f"\n=== Most Stable Surface Transitions at pH {target_pH} ===")
        
        def find_transition_point(surf1_id, surf2_id, U_start, U_end):
            """Find transition point between two surfaces in given potential range"""
            surf1 = surfs[surf1_id]
            surf2 = surfs[surf2_id]
            
            # Calculate coefficients for quadratic equation: a*U^2 + b*U + c = 0
            a = surf1['A'] - surf2['A']
            b = (surf1['B'] - surf2['B']) + (surf1['H'] - 2*surf1['O'] - surf1['e']) - (surf2['H'] - 2*surf2['O'] - surf2['e'])
            c = (surf1['E_DFT'] - surf2['E_DFT']) + const * (surf1['H'] - 2*surf1['O'] - surf2['H'] + 2*surf2['O']) * target_pH + const * (log10(surf1['conc']) - log10(surf2['conc']))
            
            # Solve quadratic equation
            if abs(a) > 1e-10:  # Quadratic case
                discriminant = b*b - 4*a*c
                if discriminant >= 0:
                    U1 = (-b + np.sqrt(discriminant)) / (2*a)
                    U2 = (-b - np.sqrt(discriminant)) / (2*a)
                    for U_sol in [U1, U2]:
                        if U_start <= U_sol <= U_end:
                            # Verify this is actually a transition point
                            E1 = dg(surf1, pH=target_pH, U=U_sol, ref_surf=surfs[ref_surf_idx])
                            E2 = dg(surf2, pH=target_pH, U=U_sol, ref_surf=surfs[ref_surf_idx])
                            if abs(E1 - E2) < 1e-6:  # Energies are equal
                                return U_sol
            else:  # Linear case
                if abs(b) > 1e-10:
                    U_sol = -c / b
                    if U_start <= U_sol <= U_end:
                        # Verify this is actually a transition point
                        E1 = dg(surf1, pH=target_pH, U=U_sol, ref_surf=surfs[ref_surf_idx])
                        E2 = dg(surf2, pH=target_pH, U=U_sol, ref_surf=surfs[ref_surf_idx])
                        if abs(E1 - E2) < 1e-6:  # Energies are equal
                            return U_sol
            return None
        
        # Sequential transition finding: start from lowest potential and find next transition
        transition_points = []
        current_surf_id = int(lowest_surfaces_pH[0])  # Most stable surface at Umin
        current_U = Umin
        
        print(f"Starting from {surfs[current_surf_id]['name']} at U = {current_U:.3f} V")
        
        while current_U < Umax:
            # Find the next transition point from current surface
            next_transition = None
            next_surf_id = None
            min_next_U = Umax + 1
            
            # Check all other surfaces for the next transition
            for other_surf_id in unique_ids_set:
                if other_surf_id != current_surf_id:
                    U_trans = find_transition_point(current_surf_id, other_surf_id, current_U + 1e-6, Umax)
                    if U_trans is not None and U_trans < min_next_U:
                        # Verify this is actually the next transition by checking stability
                        U_test = U_trans + 1e-6
                        if U_test <= Umax:
                            E_current = dg(surfs[current_surf_id], pH=target_pH, U=U_test, ref_surf=surfs[ref_surf_idx])
                            E_other = dg(surfs[other_surf_id], pH=target_pH, U=U_test, ref_surf=surfs[ref_surf_idx])
                            if E_other < E_current:  # Other surface becomes more stable
                                min_next_U = U_trans
                                next_transition = U_trans
                                next_surf_id = other_surf_id
            
            if next_transition is not None:
                transition_points.append((next_transition, current_surf_id, next_surf_id))
                print(f"  Transition at U = {next_transition:.3f} V: {surfs[current_surf_id]['name']} → {surfs[next_surf_id]['name']}")
                current_surf_id = next_surf_id
                current_U = next_transition
            else:
                break  # No more transitions found
        
        # Print transition information
        if transition_points:
            print(f"\nFound {len(transition_points)} sequential transition point(s):")
            for i, (U_trans, surf_from, surf_to) in enumerate(transition_points):
                print(f"  {i+1}. U = {U_trans:.3f} V vs SHE")
                print(f"     {surfs[surf_from]['name']} → {surfs[surf_to]['name']}")
        else:
            print("No transitions found in the potential range.")
        
        # Print stability range for each surface
        print(f"\n=== Stability Ranges at pH {target_pH} ===")
        unique_surfaces = np.unique(lowest_surfaces_pH.astype(int))
        for surf_id in unique_surfaces:
            # Find where this surface is stable
            stable_indices = np.where(lowest_surfaces_pH == surf_id)[0]
            if len(stable_indices) > 0:
                U_start = Urange[stable_indices[0]]
                U_end = Urange[stable_indices[-1]]
                print(f"  {surfs[surf_id]['name']}: {U_start:.3f} to {U_end:.3f} V vs SHE")
    
    # Calculate energy at specific pH for plotting
    all_energies = []
    
    # 1D plot: Sort unique_ids_set, unique_second_ids_set by 'e' value first, then by energy at U=0
    energies_at_U0 = [(k,surfs[k]['e'] - surfs[k]['H'] + 2*surfs[k]['O'], dg(surfs[k], pH=target_pH, U=0, ref_surf=surfs[ref_surf_idx])) for k in unique_ids_set | unique_second_ids_set]
    sorted_unique_ids = [k for k, _, _ in sorted(energies_at_U0, key=lambda x: (x[1], x[2]))]  # Sort by 'e' first, then by energy

    fig2, ax2 = plt.subplots(figsize=(args.figx, args.figy))
    ax2.axis([Umin, Umax, None, None])
    ax2.set_xlabel('Potential (V vs. RHE)')
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
        
        # Calculate how many colors go to left and right
        n_left = int(np.ceil(n_colors_1d / 2))
        n_right = n_colors_1d - n_left
        
        # Create color values avoiding the gap
        if n_left > 0:
            left = np.linspace(args.cmin_1d, left_end, n_left, endpoint=True)
        else:
            left = []
        
        if n_right > 0:
            right = np.linspace(right_start, args.cmax_1d, n_right, endpoint=True)
        else:
            right = []
        
        color_values_1d = np.concatenate([left, right]) if len(left) > 0 and len(right) > 0 else (left if len(left) > 0 else right)
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
    plt.savefig(f'{png_name}_pH{target_pH}{suffix}.png', dpi=300, bbox_inches='tight', transparent=True)
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

# ========================================
# UTILITY FUNCTIONS
# ========================================

def format_name(formula):
    """
    Format chemical formula with proper subscripts and superscripts for display.
    
    This function converts ion formulas into properly formatted strings using Unicode
    characters for subscripts (numbers) and superscripts (charges).
    
    Args:
        formula (str): Chemical formula string (e.g., "Fe2+", "SO42-", "H2O")
        
    Returns:
        str: Formatted formula with Unicode subscripts and superscripts
             (e.g., "Fe²⁺", "SO₄²⁻", "H₂O")
    
    Examples:
        >>> format_name("Fe2+")
        'Fe₂⁺'
        >>> format_name("SO42-")
        'SO₄²⁻'
        >>> format_name("H2O")
        'H₂O'
    """
    # Count charge indicators
    plus_count = formula.count('+')   # Number of positive charges
    minus_count = formula.count('-')  # Number of negative charges
    
    # Remove charge indicators to get base formula
    base_formula = formula.replace('+', '').replace('-', '')
    
    # Convert numbers to Unicode subscripts
    formatted_formula = ''
    for char in base_formula:
        if char.isdigit():
            formatted_formula += SUBSCRIPT_NUMS[char]
        else:
            formatted_formula += char
    
    # Add charge display using Unicode superscripts
    if plus_count > 0:
        if plus_count == 1:
            return formatted_formula + '⁺'
        else:
            # Multi-charge positive ion (e.g., Fe²⁺)
            superscript_num = SUPERSCRIPT_NUMS.get(str(plus_count), str(plus_count))
            return formatted_formula + superscript_num + '⁺'
    elif minus_count > 0:
        if minus_count == 1:
            return formatted_formula + '⁻'
        else:
            # Multi-charge negative ion (e.g., SO₄²⁻)
            superscript_num = SUPERSCRIPT_NUMS.get(str(minus_count), str(minus_count))
            return formatted_formula + superscript_num + '⁻'
    else:
        # Neutral species
        return formatted_formula

def format_df_for_display(df):
    """
    Format DataFrame for display with appropriate number formatting.
    
    This function formats numerical columns in a DataFrame for better readability,
    using scientific notation for concentration values and fixed decimal places
    for other numerical data.
    
    Args:
        df (pandas.DataFrame): Input DataFrame with numerical data
        
    Returns:
        pandas.DataFrame: Formatted DataFrame with string representations
    
    Notes:
        - Concentration ('conc') column uses scientific notation (e.g., 1.0e-06)
        - Other numerical columns use 2 decimal places (e.g., 1.23)
        - Non-numerical columns remain unchanged
    """
    df_copy = df.copy()
    
    # Identify numerical columns for formatting
    numeric_cols = df_copy.select_dtypes(include=[np.number]).columns
    
    for col in numeric_cols:
        if col == 'conc':
            # Use scientific notation for concentration values
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.0e}")
        else:
            # Use fixed decimal places for other numerical values
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.2f}")
    
    return df_copy

if __name__ == "__main__":
    # ========================================
    # PHYSICAL AND THERMODYNAMIC CONSTANTS
    # ========================================
    
    # Unit conversion factors
    kjmol = 96.485   # Conversion factor from kJ/mol to eV
    calmol = 23.061  # Conversion factor from cal/mol to eV

    # Fundamental physical constants
    kb = 8.617e-5    # Boltzmann constant in eV/K
    T = 298.15       # Standard temperature in K (25°C)
    const = kb * T * np.log(10)  # RT×ln(10) ≈ 0.0592 eV (for Nernst equation)
    water = 56.690/calmol  # Molar volume of water in eV units

    # ========================================
    # GAS PHASE REFERENCE ENERGIES (DFT calculated)
    # ========================================
    
    # DFT total energies for gas phase molecules (eV)
    h2 = -6.77149190   # H₂ molecule total energy
    h2o = -14.23091949 # H₂O molecule total energy

    # ========================================
    # THERMODYNAMIC CORRECTIONS FOR GAS PHASES
    # ========================================
    
    # Zero-point energy (ZPE), heat capacity (CV), and entropy (TS) corrections
    # for H₂O gas at standard conditions
    zpeh2o = 0.558  # Zero-point energy correction (eV)
    cvh2o = 0.103   # Heat capacity correction (eV)
    tsh2o = 0.675   # Entropy correction T×S (eV)

    # Corrections for H₂ gas at standard conditions
    zpeh2 = 0.268   # Zero-point energy correction (eV)
    cvh2 = 0.0905   # Heat capacity correction (eV)
    tsh2 = 0.408    # Entropy correction T×S (eV)

    # ========================================
    # GIBBS FREE ENERGIES OF GAS PHASES
    # ========================================
    
    # Calculate Gibbs free energies: G = E_DFT + ZPE + CV - TS
    gh2o = h2o + zpeh2o - tsh2o + cvh2o  # G(H₂O) gas phase
    gh2 = h2 + zpeh2 - tsh2 + cvh2       # G(H₂) gas phase

    # Derived chemical potentials for electrochemical reactions
    gh = gh2 / 2                    # μ(H) = ½G(H₂) - chemical potential of atomic H
    go = gh2o - gh2                 # μ(O) = G(H₂O) - G(H₂) - chemical potential of atomic O
    goh = gh2o - gh2 / 2           # μ(OH) = G(H₂O) - ½G(H₂) - chemical potential of OH
    gooh = 2 * gh2o - 1.5 * gh2    # μ(OOH) = 2G(H₂O) - 1.5G(H₂) - chemical potential of OOH

    # ========================================
    # ADSORBATE THERMODYNAMIC CORRECTIONS
    # ========================================
    
    # Zero-point energy, heat capacity, and entropy corrections for surface-bound species
    # OH adsorbate corrections
    zpeoh = 0.376   # Zero-point energy (eV)
    cvoh = 0.042    # Heat capacity (eV)
    tsoh = 0.066    # Entropy T×S (eV)

    # O adsorbate corrections
    zpeo = 0.064    # Zero-point energy (eV)
    cvo = 0.034     # Heat capacity (eV)
    tso = 0.060     # Entropy T×S (eV)

    # OOH adsorbate corrections
    zpeooh = 0.471  # Zero-point energy (eV)
    cvooh = 0.077   # Heat capacity (eV)
    tsooh = 0.134   # Entropy T×S (eV)

    # ========================================
    # ADSORBATE GIBBS FREE ENERGY CORRECTIONS
    # ========================================
    
    # Calculate Gibbs free energy corrections for adsorbates: ΔG = ZPE + CV - TS
    dgo = zpeo + cvo - tso        # Gibbs correction for O* adsorbate
    dgoh = zpeoh + cvoh - tsoh    # Gibbs correction for OH* adsorbate
    dgooh = zpeooh + cvooh - tsooh # Gibbs correction for OOH* adsorbate
    dgh = dgoh - dgo              # Gibbs correction for H* adsorbate (derived)

    # ========================================
    # MAIN PROGRAM EXECUTION
    # ========================================
        
    main()