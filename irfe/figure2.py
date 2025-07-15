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
adsorbates = ['*', '*OH', '*O']

results = []

def get_charge(ads, surf, path):
    atoms_path = os.path.join(root, path, 'atoms_bader_charge.json')
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

get_charge('*',   'Ir',   'slab/0_Ir')
# get_charge('*H',  'Ir',   '1_H/0_Ir/1_layer_top')
get_charge('*OH', 'Ir',   '2_OH/0_Ir/1_layer_top')
get_charge('*O',  'Ir',   '3_O/0_Ir/2_layer_hol')
get_charge('*',   'IrFe', 'slab/6_IrFe')
# get_charge('*H',  'IrFe', '1_H/6_IrFe/1_layer_top')
get_charge('*OH', 'IrFe', '2_OH/6_IrFe/1_layer_top')
get_charge('*O',  'IrFe', '3_O/6_IrFe/1_layer_top')

df = pd.DataFrame(results)
df.to_csv('bader_compact.csv', index=False)
print(df)

plt.figure(figsize=(2.5, 3))

grays = plt.cm.Greys(np.linspace(0.6, 0.2, 4))
# grays = plt.cm.Greys(np.linspace(0.9, 0.4, 4))
purples = plt.cm.Purples(np.linspace(0.9, 0.4, 4))
greens = plt.cm.Greens(np.linspace(0.9, 0.4, 4))
blues = plt.cm.Blues(np.linspace(0.9, 0.4, 4))
reds = plt.cm.Reds(np.linspace(0.8, 0.4, 3))
x = range(len(adsorbates))

# for i, layer in zip(range(4), range(4, 0, -1)):
#     charges = df[(df['layer'] == layer) & (df['surf']=='Ir')]['Ir']
#     plt.plot(x, charges, marker='o', color=grays[i], 
#     markerfacecolor='white', markeredgecolor=grays[i], 
#     linewidth=1, dashes=(3, 1), label=f'Ir L{layer}', zorder=1)

charges = df[(df['layer'] == 4) & (df['surf']=='Ir')]['Ir']
plt.plot(x, charges, marker='o', color='silver', 
         markerfacecolor='white', markeredgecolor='silver', 
         linewidth=1, dashes=(3, 1), label=f'Ir L4', zorder=1)

for i, layer in zip(range(4), range(4, 0, -1)):
    charges = df[(df['layer'] == layer) & (df['surf']=='IrFe')]['Ir']
    plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir(-Fe) L{4-layer}', zorder=2)

for i, layer in zip(range(3), range(3, 0, -1)):
    charges = df[(df['layer'] == layer) & (df['surf']=='IrFe')]['Fe']
    plt.plot(x, charges, marker='s', color=reds[i], label=f'(Ir-)Fe L{4-layer}', zorder=2)

plt.ylabel('Bader Charge (e$^{\mathrm{-}}$)')
plt.xlim(-0.5, 2.5)
plt.xticks(x, adsorbates)
plt.yticks(np.arange(-1.00, 1.26, 0.25))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5)

plot_output = os.path.join(root, 'figures', f'bader_charges.png')
plt.savefig(plot_output, bbox_inches='tight')
plt.close()

# IrFe
plt.figure(figsize=(2.5, 3))

for i, layer in zip(range(4), range(4, 0, -1)):
    charges = df[(df['layer'] == layer) & (df['surf']=='IrFe')]['Ir']
    plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir L{4-layer}')

for i, layer in zip(range(3), range(3, 0, -1)):
    charges = df[(df['layer'] == layer) & (df['surf']=='IrFe')]['Fe']
    plt.plot(x, charges, marker='s', color=reds[i], label=f'Fe L{4-layer}')

plt.xlabel('Surface Coverage')
plt.ylabel('Bader Charge (e$^{\mathrm{-}}$)')
plt.xlim(-0.5, 2.5)
plt.xticks(x, adsorbates)
plt.yticks(np.arange(-1.00, 1.26, 0.25))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5)

plot_output = os.path.join(root, 'figures', 'figure2.png')
plt.savefig(plot_output, bbox_inches='tight', dpi=300, transparent=True)
plt.close()

# Ir
plt.figure(figsize=(2.5, 3))

for i, layer in zip(range(4), range(4, 0, -1)):
    charges = df[(df['layer'] == layer) & (df['surf']=='Ir')]['Ir']
    plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir L{4-layer}')

plt.ylabel('Bader Charge (e$^{\mathrm{-}}$)')
plt.xlabel('Surface Coverage')
plt.xlim(-0.5, 2.5)
plt.xticks(x, adsorbates)
plt.yticks(np.arange(-1.00, 1.26, 0.25))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5)

plot_output = os.path.join(root, 'figures', 'figure3.png')
plt.savefig(plot_output, bbox_inches='tight', dpi=300, transparent=True)
plt.close()