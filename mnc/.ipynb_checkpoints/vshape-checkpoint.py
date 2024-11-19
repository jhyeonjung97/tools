import os
import re
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
output_dir = '/pscratch/sd/j/jiuy97/6_MNC/figures/icohp'
os.makedirs(output_dir, exist_ok=True)  # Ensure the figures directory exists
pattern = re.compile(r"\s+([\d.]+)\s+([\d.]+)")
# pattern = re.compile(r"\s+(\d+)\s+([A-Za-z]+)\s+(\d+)\s+'O'\s+(\d+)\s+[\w\(\)\-]+\s+([\d.]+)\s+([\d.]+)")

elements = {
    '3d_LS': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni'],
    '3d_HS': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni'],
    # '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
adsorption_energies = {}
icohp = {}
distance = {}
color_ranges = [plt.cm.Oranges(np.linspace(0.3, 0.9, 4)),
                plt.cm.Blues(np.linspace(0.3, 0.9, 4))]

adsorbates = ['o', 'oh']
adsorbate_symbols = ['O', 'OH']
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(elements['3d_LS'], elements['4d'], elements['5d'])]

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

def dual_scatter_by_row(elems, x1, x2, y1, y2, label1, label2, xlabel, ylabel, pngname, output_dir="."):
    plt.figure(figsize=(6.0, 4.5), dpi=300)

    x1 = np.array(x1)
    x2 = np.array(x2)
    y1 = np.array(y1)
    y2 = np.array(y2)
    
    # Scatter plots
    plt.scatter(x1, y1, color='orange', label=label1)
    plt.scatter(x2, y2, color='dodgerblue', label=label2)
    
    # Regression and R^2 calculation for the first dataset
    valid_mask1 = np.isfinite(x1) & np.isfinite(y1)
    x1_valid, y1_valid = x1[valid_mask1], y1[valid_mask1]
    z1 = np.polyfit(x1_valid, y1_valid, 1)
    p1 = np.poly1d(z1)
    y1_pred = p1(x1_valid)
    ss_res1 = np.sum((y1_valid - y1_pred) ** 2)
    ss_tot1 = np.sum((y1_valid - np.mean(y1_valid)) ** 2)
    r2_1 = 1 - (ss_res1 / ss_tot1)
    plt.plot(x1_valid, p1(x1_valid), color='orange', linestyle="--")
    equation_text1 = f"y = {z1[0]:.2f}x + {z1[1]:.2f}\n$R^2$ = {r2_1:.2f}"
    plt.text(0.05, 0.85, equation_text1, fontsize=10, color='orange', 
             ha='left', va='top', transform=plt.gca().transAxes)
    
    # Regression and R^2 calculation for the second dataset
    valid_mask2 = np.isfinite(x2) & np.isfinite(y2)
    x2_valid, y2_valid = x2[valid_mask2], y2[valid_mask2]
    z2 = np.polyfit(x2_valid, y2_valid, 1)
    p2 = np.poly1d(z2)
    y2_pred = p2(x2_valid)
    ss_res2 = np.sum((y2_valid - y2_pred) ** 2)
    ss_tot2 = np.sum((y2_valid - np.mean(y2_valid)) ** 2)
    r2_2 = 1 - (ss_res2 / ss_tot2)
    plt.plot(x2_valid, p2(x2_valid), color='dodgerblue', linestyle="--")
    equation_text2 = f"y = {z2[0]:.2f}x + {z2[1]:.2f}\n$R^2$ = {r2_2:.2f}"
    plt.text(0.05, 0.75, equation_text2, fontsize=10, color='dodgerblue', 
             ha='left', va='top', transform=plt.gca().transAxes)
    
    # Annotations
    for xi, yi, elem in zip(x1_valid, y1_valid, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='orange')
    for xi, yi, elem in zip(x2_valid, y2_valid, elems):
        plt.annotate(f'{elem}', (float(xi), float(yi)), textcoords="offset points", xytext=(0, 5), ha='center', color='dodgerblue')
    
    # Labels and legend
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    
    # Save plot
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
                       label1=f'O Adsorption ({row})', label2=f'OH Adsorption ({row})', 
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
    dual_plot_by_row(elems, x1=icohp[row]['o'], x2=icohp[row]['oh'], 
                     y1=adsorption_energies[row]['o'], y2=adsorption_energies[row]['oh'], 
                     label1=f'*O ({row})', label2=f'*OH ({row})', 
                     xlabel='-ICOHP (eV)', ylabel='Adsorption Energy (dG, eV)', pngname=f'adsorption_vs_icohp_{row}.png')
    dual_plot_by_row(elems, x1=adsorption_energies[row]['o'], x2=adsorption_energies[row]['oh'], 
                     y1=distance[row]['o'], y2=distance[row]['oh'], 
                     label1=f'*O ({row})', label2=f'*OH ({row})', 
                     xlabel='Adsorption Energy (dG, eV)', ylabel='Bond Length (Å)', pngname=f'bond_vs_adsorption_{row}.png')
    
    dual_scatter_by_row(elems, x1=distance[row]['o'], x2=distance[row]['oh'], 
                        y1=icohp[row]['o'], y2=icohp[row]['oh'], 
                        label1=f'*O ({row})', label2=f'*OH ({row})', 
                        xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', pngname=f'icohp_vs_bond_scatter_{row}.png')
    dual_scatter_by_row(elems, x1=icohp[row]['o'], x2=icohp[row]['oh'], 
                        y1=adsorption_energies[row]['o'], y2=adsorption_energies[row]['oh'], 
                        label1=f'*O ({row})', label2=f'*OH ({row})', 
                        xlabel='-ICOHP (eV)', ylabel='Adsorption Energy (dG, eV)', pngname=f'adsorption_vs_icohp_scatter_{row}.png')
    dual_scatter_by_row(elems, x1=adsorption_energies[row]['o'], x2=adsorption_energies[row]['oh'], 
                        y1=distance[row]['o'], y2=distance[row]['oh'], 
                        label1=f'*O ({row})', label2=f'*OH ({row})', 
                        xlabel='Adsorption Energy (dG, eV)', ylabel='Bond Length (Å)', pngname=f'bond_vs_adsorption_scatter_{row}.png')
    
def single_plot_by_ads(Y, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(range(len(elements[row])), Y[row][adsorbate], 
                 marker='o', color=color_ranges[i][j], label=row)
    plt.xlabel('Element')
    plt.ylabel(ylabel)
    plt.xticks(np.arange(len(indice)), indice)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
def dual_plot_by_ads(X, Y, xlabel, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        plt.plot(X[row][adsorbate], Y[row][adsorbate], marker='o', color=color_ranges[i][j], label=row)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
def dual_scatter_by_ads(X, Y, xlabel, ylabel, pngname):
    plt.figure(figsize=(6.0, 4.5), dpi=300)
    for j, row in enumerate(elements.keys()):
        xx, yy = X[row][adsorbate], Y[row][adsorbate]
        plt.scatter(xx, yy, color=color_ranges[i][j], label=row)
        valid_mask = np.isfinite(xx) & np.isfinite(yy)
        x_valid, y_valid = xx[valid_mask], yy[valid_mask]
        z = np.polyfit(x_valid, y_valid, 1)
        p = np.poly1d(z)
        y_pred = p1(x_valid)
        ss_res = np.sum((y_valid - y_pred) ** 2)
        ss_tot = np.sum((y_valid - np.mean(y_valid)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        plt.plot(x_valid, p(x_valid), color=color_ranges[i][j], linestyle="--")
        equation_text = f"y = {z[0]:.2f}x + {z[1]:.2f}\n$R^2$ = {r2:.2f}"
        plt.text(0.05, 0.85, equation_text, fontsize=10, color=color_ranges[i][j],
                 ha='left', va='top', transform=plt.gca().transAxes)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pngname), bbox_inches='tight')
    print(f"Figure saved as {pngname}")
    plt.close()
    
# Save combined plots for each adsorbate across all rows
for i in range(2): # 0, 1
    adsorbate = adsorbates[i]
    adsorbate_symbol = adsorbate_symbols[i]
    
    single_plot_by_ads(Y=adsorption_energies, ylabel=f'{adsorbate_symbol} Adsorption Energy (dG, eV)', pngname=f'adsorption_{adsorbate}.png')
    single_plot_by_ads(Y=distance, ylabel=f'M-{adsorbate_symbol} Bond Length (Å)', pngname=f'bond_{adsorbate}.png')
    single_plot_by_ads(Y=icohp, ylabel=f'-ICOHP (M-{adsorbate_symbol}, eV)', pngname=f'icohp_{adsorbate}.png')
    dual_plot_by_ads(X=distance, Y=icohp, xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', 
                     pngname=f'icohp_vs_bond_{adsorbate}.png')
    dual_plot_by_ads(X=icohp, Y=adsorption_energies, xlabel='-ICOHP (eV)', ylabel='Adsorption Energy (dG, eV)', 
                     pngname=f'adsorption_vs_icohp_{adsorbate}.png')
    dual_plot_by_ads(X=adsorption_energies, Y=distance, xlabel='Adsorption Energy (dG, eV)', ylabel='Bond Length (Å)', 
                     pngname=f'bond_vs_adsorption_{adsorbate}.png')
    dual_scatter_by_ads(X=distance, Y=icohp, xlabel='Bond Length (Å)', ylabel='-ICOHP (eV)', 
                     pngname=f'icohp_vs_bond_scatter_{adsorbate}.png')
    dual_scatter_by_ads(X=icohp, Y=adsorption_energies, xlabel='-ICOHP (eV)', ylabel='Adsorption Energy (dG, eV)', 
                     pngname=f'adsorption_vs_icohp_scatter_{adsorbate}.png')
    dual_scatter_by_ads(X=adsorption_energies, Y=distance, xlabel='Adsorption Energy (dG, eV)', ylabel='Bond Length (Å)', 
                     pngname=f'bond_vs_adsorption_scatter_{adsorbate}.png')

def data_to_tsv(data_dict, ads, filename):
    data = {
        "3d_LS": data_dict['3d_LS'][ads],
        "3d_HS": data_dict['3d_HS'][ads],
        "4d": data_dict['4d'][ads],
        "5d": data_dict['5d'][ads]
    }
    df = pd.DataFrame(data)
    df.to_csv(f"/pscratch/sd/j/jiuy97/6_MNC/figures/icohp/{filename}_{ads}.tsv", sep="\t", float_format="%.2f")
    print(f"DataFrame saved as {filename}_{ads}.tsv")

data_to_tsv(adsorption_energies, 'o', "adsorption_energies")
data_to_tsv(adsorption_energies, 'oh', "adsorption_energies")
data_to_tsv(distance, 'o', "distance")
data_to_tsv(distance, 'oh', "distance")
data_to_tsv(icohp, 'o', "icohp")
data_to_tsv(icohp, 'oh', "icohp")
