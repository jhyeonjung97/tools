import os
import json
import glob
from ase.io import read
import pandas as pd
import numpy as np

# Base directory
base_dir = '/Users/hailey/Desktop/4_IrFe3'

# Initialize results list
results = []

# Find all atoms_bader_charge.json files
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = glob.glob(pattern)

for json_file in json_files:
    # Extract path information
    path_parts = json_file.split('/')
    adsorbate = path_parts[-4]
    system_name = path_parts[-3]
    site = path_parts[-2]
    
    # Read the atoms with Bader charges
    atoms = read(json_file)
    
    # Initialize layer charges
    layer_charges = {}
    
    # Process each layer based on system_name
    if system_name == '0_Ir':
        # Original layer structure (8 atoms per layer)
        for layer in range(4):
            start_idx = layer * 8
            end_idx = (layer + 1) * 8
            layer_atoms = atoms[start_idx:end_idx]
            
            # Calculate average Bader charges for Ir and TM
            ir_charges = []
            tm_charges = []
            
            for atom in layer_atoms:
                symbol = atom.symbol
                charge = atom.charge
                
                if symbol == 'Ir':
                    ir_charges.append(charge)
                elif symbol in ['Mn', 'Fe', 'Co', 'Ni']:
                    tm_charges.append(charge)
            
            # Store average charges for this layer, rounded to 2 decimal places
            layer_charges[f'Layer{layer+1}_Ir'] = round(sum(ir_charges) / len(ir_charges), 2) if ir_charges else None
            layer_charges[f'Layer{layer+1}_TM'] = round(sum(tm_charges) / len(tm_charges), 2) if tm_charges else None
    else:
        # New layer structure
        layer_indices = [
            # Layer 1: Ir(0,1) + TM(8-13)
            {'ir': [0, 1], 'tm': list(range(8, 14))},
            # Layer 2: Ir(2,3) + TM(14-19)
            {'ir': [2, 3], 'tm': list(range(14, 20))},
            # Layer 3: Ir(4,5) + TM(20-25)
            {'ir': [4, 5], 'tm': list(range(20, 26))},
            # Layer 4: Ir(6,7) + TM(26-31)
            {'ir': [6, 7], 'tm': list(range(26, 32))}
        ]
        
        for layer_idx, indices in enumerate(layer_indices):
            # Get Ir atoms for this layer
            ir_atoms = [atoms[i] for i in indices['ir']]
            # Get TM atoms for this layer
            tm_atoms = [atoms[i] for i in indices['tm']]
            
            # Calculate Bader charges for Ir
            ir_charges = [atom.charge for atom in ir_atoms]
            
            # Calculate Bader charges for TM
            tm_charges = [atom.charge for atom in tm_atoms]
            
            # Store average charges for this layer, rounded to 2 decimal places
            layer_charges[f'Layer{layer_idx+1}_Ir'] = round(sum(ir_charges) / len(ir_charges), 2) if ir_charges else None
            layer_charges[f'Layer{layer_idx+1}_TM'] = round(sum(tm_charges) / len(tm_charges), 2) if tm_charges else None
    
    # Add to results
    result = {
        'System': system_name,
        'Site': site,
        'Adsorbate': adsorbate,
        **layer_charges
    }
    results.append(result)

# Convert to DataFrame and save as tsv
df = pd.DataFrame(results)
output_file = os.path.join(base_dir, 'bader_charges.tsv')
df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')

print(f"Results saved to {output_file}") 