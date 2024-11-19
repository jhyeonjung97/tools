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
    '3d_LS': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '3d_HS': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    # '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
adsorption_energies = {}
icohp = {}
distance = {}
color_ranges = [plt.cm.Oranges(np.linspace(0.3, 0.9, 3)),
                plt.cm.Blues(np.linspace(0.3, 0.9, 3))]

adsorbates = ['o', 'oh']
adsorbate_symbols = ['O', 'OH']
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(elements['3d'], elements['4d'], elements['5d'])]

def single_plot_by_row(elems, x1, x2, label1, label2, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.plot(elems, x1, marker='o', color='orange', label=label1)
    plt.plot(elems, x2, marker='o', color='dodgerblue', label=label2)
    plt.xlabel('Element')
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
def dual_plot_by_row(elems, x1, x2, y1, y2, label1, label2, xlabel, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.plot(x1, y1, marker='o', color='orange', label=label1)
    plt.plot(x2, y2, marker='o', color='dodgerblue', label=label2)
    for xi, yi, elem in zip(x1, y1, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='orange')
    for xi, yi, elem in zip(x2, y2, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='dodgerblue')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()

def dual_scatter_by_row(elems, x1, x2, y1, y2, label1, label2, xlabel, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    plt.scatter(x1, y1, color='orange', label=label1)
    plt.scatter(x2, y2, color='dodgerblue', label=label2)
    x = x1 + x2
    y = y1 + y2
    sorted_indices = np.argsort(x)
    x_sorted = np.array(x)[sorted_indices]
    y_sorted = np.array(y)[sorted_indices]
    z = np.polyfit(x_sorted, y_sorted, 1)
    p = np.poly1d(z)
    y_pred = p(x_sorted)
    ss_res = np.sum((y_sorted - y_pred) ** 2)
    ss_tot = np.sum((y_sorted - np.mean(y_sorted)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    plt.plot(x_sorted, p(x_sorted), color='black', linestyle="--")
    equation_text = f"y = {z[0]:.2f}x + {z[1]:.2f}\n$R^2$ = {r2:.2f}"
    plt.text(0.05, 0.15, equation_text, fontsize=10, color='black', 
             ha='left', va='top', transform=plt.gca().transAxes)
    for xi, yi, elem in zip(x1, y1, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='orange')
    for xi, yi, elem in zip(x2, y2, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='dodgerblue')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()

# Loop through each element row and calculate adsorption energies
for row, elems in elements.items():
    adsorption_energies[row] = {'o': [], 'oh': []}
    icohp[row] = {'o': [], 'oh': []}
    distance[row] = {'o': [], 'oh': []}
    for i, elem in enumerate(elems):
        numb = i + 2

        if row == '3d_LS':
            base_dir = os.path.join(root_dir, '3d', f"{numb}_{elem}", '1_LS')
        elif row == '3d_HS':
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
                if row == '3d_HS' and elem == 'Ti' and ads == 'o':
                    json_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/2_Ti/1_LS/o/final_with_calculator.json'
                elif row == '3d_HS' and elem == 'Ti' and ads == 'oh':
                    json_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/2_Ti/1_LS/oh/final_with_calculator.json'
                elif row == '3d_HS' and elem == 'V' and ads == 'o':
                    json_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/3_V/1_LS/o/final_with_calculator.json'
                else:
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
                    
                if row == '3d_HS' and elem == 'Ti' and ads == 'o':
                    icohp_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/2_Ti/1_LS/o/icohp/icohp.txt'
                elif row == '3d_HS' and elem == 'Ti' and ads == 'oh':
                    icohp_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/2_Ti/1_LS/oh/icohp/icohp.txt'
                elif row == '3d_HS' and elem == 'V' and ads == 'o':
                    icohp_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/3_V/1_LS/o/icohp/icohp.txt'
                elif row == '3d_HS' and elem == 'V' and ads == 'oh':
                    icohp_path = '/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/3_V/1_LS/oh/icohp/icohp.txt'
                else:
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
                    
    single_plot_by_row(elems, x1=adsorption_energies[row]['o'], x2=adsorption_energies[row]['oh'], 
                       label1=f'O Adsorption ({row})', label2=f'O Adsorption ({row})', 
                       ylabel='Adsorption Energy (dG, eV)', pngname=f'adsorption_{row}.png')
    single_plot_by_row(elems, x1=distance[row]['o'], x2=distance[row]['oh'], 
                       label1=f'M-O Bond ({row})', label2=f'M-OH Bond ({row})', 
                       ylabel='Bond Length (Å)', pngname=f'bond_{row}.png')
    single_plot_by_row(elems, x1=icohp[row]['o'], x2=icohp[row]['oh'], 
                       label1=f'*O ({row})', label2=f'*OH ({row})', 
                       ylabel='-ICOHP (eV)', pngname=f'icohp_{row}.png')
    dual_plot_by_row(elems, x1=distance[row]['o'], x2=distance[row]['oh'], 
                     y1=icohp[row]['o'], y2=icohp[row]['oh'], 
                     label1=f'*O ({row})', label2=f'*OH ({row})', 
                     xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', pngname=f'icohp_vs_bond_{row}.png')
    dual_scatter_by_row(elems, x1=distance[row]['o'], x2=distance[row]['oh'], 
                        y1=icohp[row]['o'], y2=icohp[row]['oh'], 
                        label1=f'*O ({row})', label2=f'*OH ({row})', 
                        xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', pngname=f'icohp_vs_bond_scatter_{row}.png')

def single_plot_by_ads(Y, label, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(range(len(elements[row])), Y[row][adsorbate], 
                 marker='o', color=color_ranges[i][j], label=label)
    plt.xlabel('Element')
    plt.ylabel(ylabel)
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
def dual_plot_by_ads(X, Y, label, xlabel, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(X[row][adsorbate], Y[row][adsorbate], 
                 marker='o', color=color_ranges[i][j], label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
# Save combined plots for each adsorbate across all rows
for i in range(2): # 0, 1
    adsorbate = adsorbates[i]
    adsorbate_symbol = adsorbate_symbols[i]
    
    single_plot_by_ads(Y=adsorption_energies, label=f'{adsorbate_symbol} Adsorption ({row})', 
                       ylabel='Adsorption Energy (dG, eV)', pngname=f'adsorption_{adsorbate}.png')
    single_plot_by_ads(Y=distance, label=f'M-{adsorbate_symbol} Bond ({row})', 
                       ylabel='Bond Length (Å)', pngname=f'bond_{adsorbate}.png')
    single_plot_by_ads(Y=icohp, label=f'*{adsorbate_symbol} ({row})', 
                       ylabel='-ICOHP (eV)', pngname=f'icohp_{adsorbate}.png')
    dual_plot_by_ads(X=distance, Y=icohp, label=f'*{adsorbate_symbol} ({row})', 
                     xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', pngname=f'icohp_vs_bond_{adsorbate}.png')