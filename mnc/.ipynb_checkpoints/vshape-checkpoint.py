import os
import re
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
pattern = re.compile(r"\s+([\d.]+)\s+([\d.]+)")
# pattern = re.compile(r"\s+(\d+)\s+([A-Za-z]+)\s+(\d+)\s+'O'\s+(\d+)\s+[\w\(\)\-]+\s+([\d.]+)\s+([\d.]+)")

elements = {
    # '3d(LS)': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    # '3d(HS)': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
adsorption_energies = {}
icohp = {}
distance = {}
color_ranges = [plt.cm.Oranges(np.linspace(0.3, 0.9, 3)),
                plt.cm.Blues(np.linspace(0.3, 0.9, 3))]

# Loop through each element row and calculate adsorption energies
for row, elems in elements.items():
    adsorption_energies[row] = {'o': [], 'oh': []}
    icohp[row] = {'o': [], 'oh': []}
    distance[row] = {'o': [], 'oh': []}
    for i, elem in enumerate(elems):
        numb = i + 2

        if row == '3d(LS)':
            base_dir = os.path.join(root_dir, '3d', f"{numb}_{elem}", '1_LS')
        elif row == '3d(HS)':
            base_dir = os.path.join(root_dir, '3d', f"{numb}_{elem}", '2_HS')
            if not os.path.exists(base_dir):
                base_dir = os.path.join(root_dir, '3d', f"{numb}_{elem}", '3_HS')
        else:
            base_dir = os.path.join(root_dir, row, f"{numb}_{elem}", 'most_stable')

        relaxed_dir = os.path.join(base_dir, 'relaxed', 'final_with_calculator.json')
        if os.path.exists(relaxed_dir):
            relaxed_atoms = read(relaxed_dir)
            e_clean = relaxed_atoms.get_potential_energy()

            for ads in ['o', 'oh']:
                json_path = os.path.join(base_dir, ads, 'final_with_calculator.json')
                if os.path.exists(json_path):        
                    atoms = read(json_path)
                    e_ads = atoms.get_potential_energy()
                    if ads == 'o':
                        dg_ads = (e_ads + dgo) - e_clean - go
                    elif ads == 'oh':
                        dg_ads = (e_ads + dgoh) - e_clean - goh
                    adsorption_energies[row][ads].append(dg_ads)
                else:
                    adsorption_energies[row][ads].append(np.nan)
                    
                icohp_path = os.path.join(base_dir, ads, 'icohp/icohp.txt')
                if os.path.exists(icohp_path):        
                    with open(icohp_path, 'r') as file:
                        content = file.read()
                    matches = pattern.findall(content)
                    for match in matches:
                        icohp[row][ads].append(float(match[0]))
                        distance[row][ads].append(float(match[1]))
                else:
                    icohp[row][ads].append(np.nan)  # Handle missing data
                    distance[row][ads].append(np.nan)  # Handle missing data

    # Save individual plots for each row
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.plot(elems, adsorption_energies[row]['o'], marker='o', color='orange', label=f'O Adsorption ({row})')
    plt.plot(elems, adsorption_energies[row]['oh'], marker='o', color='dodgerblue', label=f'OH Adsorption ({row})')
    plt.xlabel('Element')
    plt.ylabel('Adsorption Energy (dG, eV)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'adsorption_{row}.png'), bbox_inches='tight')
    print(f"Figure saved as adsorption_{row}.png")
    plt.close()
    
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.plot(elems, distance[row]['o'], marker='o', color='orange', label=f'M-O Bond ({row})')
    plt.plot(elems, distance[row]['oh'], marker='o', color='dodgerblue', label=f'M-OH Bond ({row})')
    plt.xlabel('Element')
    plt.ylabel('Bond Length (Å)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'bond_{row}.png'), bbox_inches='tight')
    print(f"Figure saved as bond_{row}.png")
    plt.close()
    
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.plot(elems, icohp[row]['o'], marker='o', color='orange', label=f'*O ({row})')
    plt.plot(elems, icohp[row]['oh'], marker='o', color='dodgerblue', label=f'*OH ({row})')
    plt.xlabel('Element')
    plt.ylabel('-ICOHP (eV)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'icohp_{row}.png'), bbox_inches='tight')
    print(f"Figure saved as icohp_{row}.png")
    plt.close()

    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.scatter(distance[row]['o'], icohp[row]['o'], color='orange', label=f'*O ({row})')
    for xi, yi, elem in zip(distance[row]['o'], icohp[row]['o'], elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='orange')
    plt.scatter(distance[row]['oh'], icohp[row]['oh'], color='dodgerblue', label=f'*OH ({row})')
    for xi, yi, elem in zip(distance[row]['oh'], icohp[row]['oh'], elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='dodgerblue')
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('-ICOHP (eV)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'icohp_vs_bond_{row}.png'), bbox_inches='tight')
    print(f"Figure saved as icohp_vs_bond_{row}.png")
    plt.close()
    
# Save combined plots for each adsorbate across all rows
adsorbates = ['o', 'oh']
adsorbate_symbols = ['O', 'OH']
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(elements['3d'], elements['4d'], elements['5d'])]

for i in range(2): # 0, 1
    adsorbate = adsorbates[i]
    adsorbate_symbol = adsorbate_symbols[i]
    
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(range(len(elements[row])), adsorption_energies[row][adsorbate], 
                 marker='o', color=color_ranges[i][j],
                 label=f'{adsorbate_symbol} Adsorption ({row})')
    plt.xlabel('Element')
    plt.ylabel('Adsorption Energy (dG, eV)')
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'adsorption_{adsorbate}.png'), bbox_inches='tight')
    print(f"Figure saved as adsorption_{adsorbate}.png")
    plt.close()
    
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(range(len(elements[row])), distance[row][adsorbate], 
                 marker='o', color=color_ranges[i][j],
                 label=f'{adsorbate_symbol} Adsorption ({row})')
    plt.xlabel('Element')
    plt.ylabel('Bond Length (Å)')
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'bond_{adsorbate}.png'), bbox_inches='tight')
    print(f"Figure saved as bond_{adsorbate}.png")
    plt.close()

    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(range(len(elements[row])), icohp[row][adsorbate], 
                 marker='o', color=color_ranges[i][j],
                 label=f'*{adsorbate_symbol} ({row})')
    plt.xlabel('Element')
    plt.ylabel('-ICOHP (eV)')
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'icohp_{adsorbate}.png'), bbox_inches='tight')
    print(f"Figure saved as icohp_{adsorbate}.png")
    plt.close()

    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.scatter(range(len(elements[row])), distance[row][adsorbate], 
                 color=color_ranges[i][j], label=f'*{adsorbate_symbol} ({row})')
        for xi, yi, elem in zip(range(len(elements[row])), distance[row][adsorbate], elements[row]):
            plt.annotate(f'{elem}', (float(xi), float(yi)), 
                         textcoords="offset points", xytext=(0, 5), ha='center', color=color_ranges[i][j])
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('-ICOHP (eV)')
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'icohp_vs_bond_{adsorbate}.png'), bbox_inches='tight')
    print(f"Figure saved as icohp_vs_bond_{adsorbate}.png")
    plt.close()
