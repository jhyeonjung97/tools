import os
import json
import glob
from ase.io import read
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.geometry import get_distances

# Base directory
base_dir = '/Users/hailey/Desktop/4_IrFe3'

# Initialize results list
results = []

# Find all atoms_bader_charge.json files and sort them
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = sorted(glob.glob(pattern))

# Find slab atoms_bader_charge.json files
slab_pattern = os.path.join(base_dir, 'slab/*_*/atoms_bader_charge.json')
slab_json_files = sorted(glob.glob(slab_pattern))

def get_layer_charges_and_distances(atoms):
    # Get z-coordinates and sort atoms
    z_coords = atoms.get_positions()[:, 2]
    sorted_indices = np.argsort(z_coords)
    sorted_atoms = atoms[sorted_indices]
    
    # Group atoms by z-coordinate to identify layers
    z_coords_sorted = z_coords[sorted_indices]
    z_diff = np.diff(z_coords_sorted)
    layer_boundaries = np.where(z_diff > 0.5)[0] + 1  # Threshold for layer separation
    
    # Split atoms into layers
    layer_indices = np.split(np.arange(len(atoms)), layer_boundaries)
    
    # Initialize layer charges and distances
    layer_charges = {}
    layer_distances = {}
    
    # Process each layer
    for layer_idx, indices in enumerate(layer_indices):
        layer_atoms = sorted_atoms[indices]
        
        # Calculate average Bader charges for Ir and TM
        ir_charges = []
        tm_charges = []
        ir_indices = []
        tm_indices = []
        
        for i, atom in enumerate(layer_atoms):
            symbol = atom.symbol
            charge = atom.charge
            
            if symbol == 'Ir':
                ir_charges.append(charge)
                ir_indices.append(i)
            elif symbol in ['Mn', 'Fe', 'Co', 'Ni']:
                tm_charges.append(charge)
                tm_indices.append(i)
        
        # Store average charges for this layer
        layer_charges[f'L{layer_idx+1}_Ir'] = round(sum(ir_charges) / len(ir_charges), 2) if ir_charges else None
        layer_charges[f'L{layer_idx+1}_TM'] = round(sum(tm_charges) / len(tm_charges), 2) if tm_charges else None
        
        # Calculate distances within the layer
        if ir_indices and tm_indices:
            # Calculate Ir-TM distances
            ir_positions = layer_atoms[ir_indices].get_positions()
            tm_positions = layer_atoms[tm_indices].get_positions()
            distances = get_distances(ir_positions, tm_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Find all distances within threshold
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_idx+1}_Ir-TM'] = round(np.mean(valid_distances), 2)
        
        if len(ir_indices) > 1:
            # Calculate Ir-Ir distances within layer
            ir_positions = layer_atoms[ir_indices].get_positions()
            distances = get_distances(ir_positions, ir_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Get non-zero distances within threshold
            distances = distances[distances > 0]
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_idx+1}_Ir-Ir'] = round(np.mean(valid_distances), 2)
        
        if len(tm_indices) > 1:
            # Calculate TM-TM distances within layer
            tm_positions = layer_atoms[tm_indices].get_positions()
            distances = get_distances(tm_positions, tm_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Get non-zero distances within threshold
            distances = distances[distances > 0]
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_idx+1}_TM-TM'] = round(np.mean(valid_distances), 2)
        
        # Calculate distances between layers
        if layer_idx < len(layer_indices) - 1:
            next_layer_atoms = sorted_atoms[layer_indices[layer_idx + 1]]
            next_ir_indices = [i for i, atom in enumerate(next_layer_atoms) if atom.symbol == 'Ir']
            next_tm_indices = [i for i, atom in enumerate(next_layer_atoms) if atom.symbol in ['Mn', 'Fe', 'Co', 'Ni']]
            
            # Calculate Ir-TM distances between layers
            if ir_indices and next_tm_indices:
                ir_positions = layer_atoms[ir_indices].get_positions()
                next_tm_positions = next_layer_atoms[next_tm_indices].get_positions()
                distances = get_distances(ir_positions, next_tm_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_idx+1}-L{layer_idx+2}_Ir-TM'] = round(np.mean(valid_distances), 2)
            
            if tm_indices and next_ir_indices:
                tm_positions = layer_atoms[tm_indices].get_positions()
                next_ir_positions = next_layer_atoms[next_ir_indices].get_positions()
                distances = get_distances(tm_positions, next_ir_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_idx+1}-L{layer_idx+2}_TM-Ir'] = round(np.mean(valid_distances), 2)
            
            if ir_indices and next_ir_indices:
                # Calculate Ir-Ir distances between layers
                ir_positions = layer_atoms[ir_indices].get_positions()
                next_ir_positions = next_layer_atoms[next_ir_indices].get_positions()
                distances = get_distances(ir_positions, next_ir_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_idx+1}-L{layer_idx+2}_Ir-Ir'] = round(np.mean(valid_distances), 2)
            
            if tm_indices and next_tm_indices:
                # Calculate TM-TM distances between layers
                tm_positions = layer_atoms[tm_indices].get_positions()
                next_tm_positions = next_layer_atoms[next_tm_indices].get_positions()
                distances = get_distances(tm_positions, next_tm_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_idx+1}-L{layer_idx+2}_TM-TM'] = round(np.mean(valid_distances), 2)
    
    return layer_charges, layer_distances

# Process regular files
for json_file in json_files:
    # Extract path information
    path_parts = json_file.split('/')
    adsorbate = path_parts[-4]
    system_name = path_parts[-3]
    site = path_parts[-2]
    
    # Read the atoms with Bader charges
    atoms = read(json_file)
    
    # Get layer charges and distances based on z-coordinates
    layer_charges, layer_distances = get_layer_charges_and_distances(atoms)
    
    # Add to results
    result = {
        'Adsorbate': adsorbate,
        'System': system_name,
        'Site': site,
        **layer_charges,
        **layer_distances
    }
    results.append(result)

# Process slab files
for json_file in slab_json_files:
    # Extract path information
    path_parts = json_file.split('/')
    system_name = path_parts[-2]  # e.g., IrFe3_111
    
    # Read the atoms with Bader charges
    atoms = read(json_file)
    
    # Get layer charges and distances based on z-coordinates
    layer_charges, layer_distances = get_layer_charges_and_distances(atoms)
    
    # Add to results
    result = {
        'Adsorbate': 'slab',
        'System': system_name,
        'Site': 'slab',
        **layer_charges,
        **layer_distances
    }
    results.append(result)

# Convert to DataFrame
df = pd.DataFrame(results)

# Save as tsv with 2 decimal places
tsv_output = os.path.join(base_dir, 'bader_charges_and_distances.tsv')
df.to_csv(tsv_output, sep='\t', index=False, float_format='%.2f')

# Save as csv with 4 decimal places
csv_output = os.path.join(base_dir, 'bader_charges_and_distances.csv')
df.to_csv(csv_output, index=False, float_format='%.4f')

print(f"Results saved to {tsv_output} and {csv_output}")

# Create plots for each adsorbate and site combination
sns.set_theme(style="ticks")

# Define colormaps for different distance types
ir_tm_cmap = plt.cm.Greens(np.linspace(0.3, 0.9, 4))  # Green gradient for Ir-TM
ir_ir_cmap = plt.cm.Blues(np.linspace(0.3, 0.9, 4))   # Blue gradient for Ir-Ir
tm_tm_cmap = plt.cm.Reds(np.linspace(0.3, 0.9, 4))    # Red gradient for TM-TM

# Get unique adsorbates and sites
adsorbates = df['Adsorbate'].unique()
sites = df['Site'].unique()

# Define system names for x-axis
system_names = ['Ir', 'IrMn', 'IrFe', 'IrCo', 'IrNi']

# Create a figure for each adsorbate-site combination
for adsorbate in adsorbates:
    for site in sites:
        # Skip slab data in the main loop
        if adsorbate == 'slab' and site == 'slab':
            continue
            
        # Filter data for current adsorbate and site
        subset = df[(df['Adsorbate'] == adsorbate) & (df['Site'] == site)]
        if len(subset) == 0:
            continue
            
        # Create figure for Bader charges
        plt.figure(figsize=(4, 3))
        
        # Prepare data for plotting
        systems = subset['System']
        x = range(len(systems))
        
        # Plot Ir charges
        for layer in range(1, 5):
            charges = subset[f'L{layer}_Ir']
            plt.plot(x, charges, marker='o', label=f'L{layer} Ir', color=ir_ir_cmap[layer-1])
        
        # Plot TM charges
        for layer in range(1, 5):
            charges = subset[f'L{layer}_TM']
            plt.plot(x, charges, marker='s', label=f'L{layer} TM', color=tm_tm_cmap[layer-1])
        
        # Customize the plot
        plt.ylabel('Bader Charge')
        plt.xticks(x, system_names)
        plt.yticks(np.arange(-1.25, 1.01, 0.25))
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        
        # Save the plot
        plot_output = os.path.join(base_dir, 'figures', f'bader_charges_{adsorbate}_{site}.png')
        plt.savefig(plot_output, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to {plot_output}")
        
        # Create figure for distances
        plt.figure(figsize=(4, 3))
        
        # Plot intra-layer distances
        for layer in range(1, 5):
            # Ir-TM distances
            if f'L{layer}_Ir-TM' in subset.columns:
                distances = subset[f'L{layer}_Ir-TM']
                plt.plot(x, distances, marker='o', label=f'L{layer} Ir-TM', color=ir_tm_cmap[layer-1])
            
            # Ir-Ir distances
            if f'L{layer}_Ir-Ir' in subset.columns:
                distances = subset[f'L{layer}_Ir-Ir']
                plt.plot(x, distances, marker='o', label=f'L{layer} Ir-Ir', color=ir_ir_cmap[layer-1])
            
            # TM-TM distances
            if f'L{layer}_TM-TM' in subset.columns:
                distances = subset[f'L{layer}_TM-TM']
                plt.plot(x, distances, marker='o', label=f'L{layer} TM-TM', color=tm_tm_cmap[layer-1])
        
        # Plot inter-layer distances
        for layer in range(1, 4):
            # Ir-TM distances
            if f'L{layer}-L{layer+1}_Ir-TM' in subset.columns:
                distances = subset[f'L{layer}-L{layer+1}_Ir-TM']
                plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} Ir-TM', color=ir_tm_cmap[layer-1])
            if f'L{layer}-L{layer+1}_TM-Ir' in subset.columns:
                distances = subset[f'L{layer}-L{layer+1}_TM-Ir']
                plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} TM-Ir', color=ir_tm_cmap[layer-1])
            
            # Ir-Ir distances
            if f'L{layer}-L{layer+1}_Ir-Ir' in subset.columns:
                distances = subset[f'L{layer}-L{layer+1}_Ir-Ir']
                plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} Ir-Ir', color=ir_ir_cmap[layer-1])
            
            # TM-TM distances
            if f'L{layer}-L{layer+1}_TM-TM' in subset.columns:
                distances = subset[f'L{layer}-L{layer+1}_TM-TM']
                plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} TM-TM', color=tm_tm_cmap[layer-1])
        
        # Customize the plot
        plt.ylabel('Distance (Å)')
        plt.xticks(x, system_names)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', ncol=2)
        
        # Save the plot
        plot_output = os.path.join(base_dir, 'figures', f'distances_{adsorbate}_{site}.png')
        plt.savefig(plot_output, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to {plot_output}")

# Create separate plots for slab data
slab_data = df[(df['Adsorbate'] == 'slab') & (df['Site'] == 'slab')]
if len(slab_data) > 0:
    # Create figure for Bader charges
    plt.figure(figsize=(4, 3))
    
    # Prepare data for plotting
    systems = slab_data['System']
    x = range(len(systems))
    
    # Plot Ir charges
    for layer in range(1, 5):
        charges = slab_data[f'L{layer}_Ir']
        plt.plot(x, charges, marker='o', label=f'L{layer} Ir', color=ir_ir_cmap[layer-1])
    
    # Plot TM charges
    for layer in range(1, 5):
        charges = slab_data[f'L{layer}_TM']
        plt.plot(x, charges, marker='s', label=f'L{layer} TM', color=tm_tm_cmap[layer-1])
    
    # Customize the plot
    plt.ylabel('Bader Charge')
    plt.xticks(x, system_names)
    plt.yticks(np.arange(-1.25, 1.01, 0.25))
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    
    # Save the plot
    slab_plot_output = os.path.join(base_dir, 'figures', 'bader_charges_slab.png')
    plt.savefig(slab_plot_output, bbox_inches='tight')
    plt.close()
    
    print(f"Slab plot saved to {slab_plot_output}")
    
    # Create figure for distances
    plt.figure(figsize=(4, 3))
    
    # Plot intra-layer distances
    for layer in range(1, 5):
        # Ir-TM distances
        if f'L{layer}_Ir-TM' in slab_data.columns:
            distances = slab_data[f'L{layer}_Ir-TM']
            plt.plot(x, distances, marker='o', label=f'L{layer} Ir-TM', color=ir_tm_cmap[layer-1])
        
        # Ir-Ir distances
        if f'L{layer}_Ir-Ir' in slab_data.columns:
            distances = slab_data[f'L{layer}_Ir-Ir']
            plt.plot(x, distances, marker='o', label=f'L{layer} Ir-Ir', color=ir_ir_cmap[layer-1])
        
        # TM-TM distances
        if f'L{layer}_TM-TM' in slab_data.columns:
            distances = slab_data[f'L{layer}_TM-TM']
            plt.plot(x, distances, marker='o', label=f'L{layer} TM-TM', color=tm_tm_cmap[layer-1])
    
    # Plot inter-layer distances
    for layer in range(1, 4):
        # Ir-TM distances
        if f'L{layer}-L{layer+1}_Ir-TM' in slab_data.columns:
            distances = slab_data[f'L{layer}-L{layer+1}_Ir-TM']
            plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} Ir-TM', color=ir_tm_cmap[layer-1])
        if f'L{layer}-L{layer+1}_TM-Ir' in slab_data.columns:
            distances = slab_data[f'L{layer}-L{layer+1}_TM-Ir']
            plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} TM-Ir', color=ir_tm_cmap[layer-1])
        
        # Ir-Ir distances
        if f'L{layer}-L{layer+1}_Ir-Ir' in slab_data.columns:
            distances = slab_data[f'L{layer}-L{layer+1}_Ir-Ir']
            plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} Ir-Ir', color=ir_ir_cmap[layer-1])
        
        # TM-TM distances
        if f'L{layer}-L{layer+1}_TM-TM' in slab_data.columns:
            distances = slab_data[f'L{layer}-L{layer+1}_TM-TM']
            plt.plot(x, distances, marker='^', label=f'L{layer}-L{layer+1} TM-TM', color=tm_tm_cmap[layer-1])
    
    # Customize the plot
    plt.ylabel('Distance (Å)')
    plt.xticks(x, system_names)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', ncol=2)
    
    # Save the plot
    slab_plot_output = os.path.join(base_dir, 'figures', 'distances_slab.png')
    plt.savefig(slab_plot_output, bbox_inches='tight')
    plt.close()
    
    print(f"Slab distances plot saved to {slab_plot_output}") 