import os
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np

# Gas phase energies and vibrational corrections
h2 = -6.77149190
h2o = -14.23091949
zpeh2o = 0.560
cvh2o = 0.103
tsh2o = 0.675
zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + cvh2o - tsh2o + zpeh2o
gh2 = h2 + cvh2 - tsh2 + zpeh2
go = gh2o - gh2
goh = gh2o - gh2 / 2

# Adsorbate vibrational corrections
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066
zpeo = 0.064
cvo = 0.034
tso = 0.060

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh

# Define base directory, elements, and adsorbates
root_dir = '/pscratch/sd/j/jiuy97/6_MNC/0_clean'
output_dir = '/pscratch/sd/j/jiuy97/6_MNC/figures'
os.makedirs(output_dir, exist_ok=True)  # Ensure the figures directory exists

elements = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
adsorption_energies = {}
color_ranges = [plt.cm.Oranges(np.linspace(0.3, 0.9, 3)),
                plt.cm.Blues(np.linspace(0.3, 0.9, 3))]

# Loop through each element row and calculate adsorption energies
for row, elems in elements.items():
    adsorption_energies[row] = {'o': [], 'oh': []}
    for i, elem in enumerate(elems):
        numb = i + 2
        base_dir = os.path.join(root_dir, row, f"{numb}_{elem}", 'most_stable')

        o_dir = os.path.join(base_dir, 'o', 'final_with_calculator.json')
        oh_dir = os.path.join(base_dir, 'oh', 'final_with_calculator.json')
        relaxed_dir = os.path.join(base_dir, 'relaxed', 'final_with_calculator.json')

        # Check if clean surface file exists
        if os.path.exists(relaxed_dir):
            relaxed_atoms = read(relaxed_dir)
            e_clean = relaxed_atoms.get_potential_energy()
                
            # Calculate O adsorption energy if file exists
            if os.path.exists(o_dir):        
                o_atoms = read(o_dir)
                e_o = o_atoms.get_potential_energy()
                dg_o = (e_o + dgo) - e_clean - go
                adsorption_energies[row]['o'].append(dg_o)
            else:
                adsorption_energies[row]['o'].append(None)  # Handle missing data
                
            # Calculate OH adsorption energy if file exists
            if os.path.exists(oh_dir):        
                oh_atoms = read(oh_dir)
                e_oh = oh_atoms.get_potential_energy()
                dg_oh = (e_oh + dgoh) - e_clean - goh
                adsorption_energies[row]['oh'].append(dg_oh)
            else:
                adsorption_energies[row]['oh'].append(None)  # Handle missing data

    # Save individual plots for each row
    plt.figure(figsize=(4, 3), dpi=300)
    plt.plot(elems, adsorption_energies[row]['o'], marker='o', color='orange', label=f'O Adsorption ({row})')
    plt.plot(elems, adsorption_energies[row]['oh'], marker='o', color='blue', label=f'OH Adsorption ({row})')
    plt.xlabel('Element')
    plt.ylabel('Adsorption Energy (dG, eV)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'adsorption_{row}.png'), bbox_inches='tight')
    print(f"Figure saved as adsorption_{row}.png")

# Save combined plots for each adsorbate across all rows
small_adsorbates = ['o', 'oh']
large_adsorbates = ['O', 'OH']
for i in range(2): # 0, 1
    plt.figure(figsize=(4, 3), dpi=300)
    small_adsorbate = small_adsorbates[i]
    large_adsorbate = large_adsorbates[i]
    for j, row in enumerate(elements.keys()):
        plt.plot(elements[row], adsorption_energies[row][small_adsorbate], 
                 marker='o', color=color_ranges[i][j],
                 label=f'{large_adsorbate} Adsorption ({row})')
    plt.xlabel('Element')
    plt.ylabel('Adsorption Energy (dG, eV)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'adsorption_{small_adsorbate}.png'), bbox_inches='tight')
    print(f"Figure saved as adsorption_{small_adsorbate}.png")
