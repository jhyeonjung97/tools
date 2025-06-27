import os
import json
import glob
from ase.io import read
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.geometry import get_distances
from collections import defaultdict

# Base directory
root = '/Users/hailey/Desktop/4_IrFe3'

surfaces = ['Ir(Ir)', 'Ir(Fe)', 'Ir-Fe(Ir)', 'Ir-Fe(Fe)']
adsorbates = ['*H', '*', '*OH', '*O']

results = []

def get_charge(ads, surf, path):
    atoms_path = os.path.join(root, path, 'atoms_bader_charge.json')
    print(atoms_path)
    atoms = read(atoms_path)

    z_coords = atoms.get_positions()[:, 2]
    sorted_indices = np.argsort(z_coords)
    sorted_atoms = atoms[sorted_indices]
    
    num_atoms = 32
    num_layers = 4
    layer_size = 8

    layer_indices = []
    for i in range(num_layers):
        start_idx = i * layer_size
        end_idx = min((i + 1) * layer_size, num_atoms)
        layer_indices.append(np.arange(start_idx, end_idx))
    
    print(layer_indices)
    
    for i, indices in enumerate(layer_indices):
        layer_num = i + 1
        layer_atoms = sorted_atoms[indices]
        
        ir_charges = []
        tm_charges = []
        ir_indices = []
        tm_indices = []
        
        for i, atom in enumerate(layer_atoms):
            symbol = atom.symbol
            charge = -atom.charge
            
            if symbol == 'Ir':
                ir_charges.append(charge)
                ir_indices.append(i)
            elif symbol in ['Mn', 'Fe', 'Co', 'Ni']:
                tm_charges.append(charge)
                tm_indices.append(i)
        
        # Store average charges for this layer
        layer_charges_Ir = round(sum(ir_charges) / len(ir_charges), 2) if ir_charges else None
        layer_charges_Fe = round(sum(tm_charges) / len(tm_charges), 2) if tm_charges else None

        results.append({
            'ads': ads,
            'surf': surf,
            'layer': layer_num,
            'Ir': layer_charges_Ir,
            'Fe': layer_charges_Fe
        })

get_charge('*H',  'Ir',   '1_H/0_Ir/1_layer_top')
get_charge('*',   'Ir',   'slab/0_Ir')
get_charge('*OH', 'Ir',   '2_OH/0_Ir/1_layer_top')
get_charge('*O',  'Ir',   '3_O/0_Ir/2_layer_hol')
get_charge('*H',  'IrFe', '1_H/6_IrFe/1_layer_top')
get_charge('*',   'IrFe', 'slab/6_IrFe')
get_charge('*OH', 'IrFe', '2_OH/6_IrFe/1_layer_top')
get_charge('*O',  'IrFe', '3_O/6_IrFe/1_layer_top')

df = pd.DataFrame(results)
print(df)
df.to_csv('bader_compact.csv', index=False)