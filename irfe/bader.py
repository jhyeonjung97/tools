import os
import json
import glob
from ase.io import read
import pandas as pd

# Base directory
# base_dir = '/home/hyeonjung/scratch/4_IrFe3'
base_dir = '/Users/hailey/Desktop/4_IrFe3'

# ZVAL values from POTCAR
zval_dict = {
    'Ir': 9,  # Ir: [Xe] 4f14 5d7 6s2
    'Mn': 7,  # Mn: [Ar] 3d5 4s2
    'Fe': 8,  # Fe: [Ar] 3d6 4s2
    'Co': 9,  # Co: [Ar] 3d7 4s2
    'Ni': 10, # Ni: [Ar] 3d8 4s2
    'O': 6,   # O: [He] 2s2 2p4
    'H': 1    # H: 1s1
}

# Initialize results list
results = []

# Find all atoms_bader_charge.json files
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = glob.glob(pattern)

for json_file in json_files:
    # Extract path information
    path_parts = json_file.split('/')
    system_name = path_parts[-4]  # e.g., IrFe3_111
    adsorbate = path_parts[-3]    # e.g., O2
    site = path_parts[-2]         # e.g., fcc
    calculation = path_parts[-3]  # e.g., O2_fcc
    
    # Read the atoms with Bader charges
    atoms = read(json_file)
    
    # Initialize layer charges
    layer_charges = {}
    
    # Process each layer (4 layers total)
    for layer in range(4):
        # Get indices for this layer (8 atoms per layer)
        start_idx = layer * 8
        end_idx = (layer + 1) * 8
        
        # Get atoms in this layer
        layer_atoms = atoms[start_idx:end_idx]
        
        # Calculate average Bader charges for Ir and TM
        ir_charges = []
        tm_charges = []
        
        for atom in layer_atoms:
            symbol = atom.symbol
            charge = atom.charge
            zval = zval_dict.get(symbol, 0)
            bader_charge = zval - charge
            
            if symbol == 'Ir':
                ir_charges.append(bader_charge)
            elif symbol in ['Mn', 'Fe', 'Co', 'Ni']:
                tm_charges.append(bader_charge)
        
        # Store average charges for this layer
        layer_charges[f'Layer{layer+1}_Ir'] = sum(ir_charges) / len(ir_charges) if ir_charges else None
        layer_charges[f'Layer{layer+1}_TM'] = sum(tm_charges) / len(tm_charges) if tm_charges else None
    
    # Add to results
    result = {
        'System': system_name,
        'Adsorbate': adsorbate,
        'Site': site,
        'Calculation': calculation,
        **layer_charges
    }
    results.append(result)

# Convert to DataFrame and save as tsv
df = pd.DataFrame(results)
output_file = os.path.join(base_dir, 'bader_charges.tsv')
df.to_csv(output_file, sep='\t', index=False)

print(f"Results saved to {output_file}") 