import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

png_filename = "energy_norm_formation.png"
tsv_filename = "energy_norm_formation.tsv"

print(f"\033[92m{os.getcwd()}\033[0m")
if '1_Tetrahedral_WZ' in os.getcwd():
    marker = '>'; color = '#d62728'; coordination = 'WZ'; n = 1
elif '2_Tetrahedral_ZB' in os.getcwd():
    marker = '<'; color = '#ff7f0e'; coordination = 'ZB'; n = 1
elif '3_SquarePlanar_TN' in os.getcwd():
    marker = 'o'; color = '#ffd70e'; coordination = 'TN'; n = 1
elif '4_SquarePlanar_PD' in os.getcwd():
    marker = 's'; color = '#2ca02c'; coordination = 'PD'; n = 1
elif '5_SquarePlanar_NB' in os.getcwd():
    marker = 'p'; color = '#17becf'; coordination = 'NB'; n = 1
elif '6_Octahedral_RS' in os.getcwd():
    marker = 'd'; color = '#9467bd'; coordination = 'RS'; n = 1
elif '7_Pyramidal_LT' in os.getcwd():
    marker = 'h'; color = '#8c564b'; coordination = 'LT'; n = 1
elif '8_Tetrahedral_AQ' in os.getcwd():
    marker = '^'; color = '#e377c2'; coordination = 'AQ'; n = 2
elif '9_SquarePlanar_AU' in os.getcwd():
    marker = 'v'; color = '#7f7f7f'; coordination = 'AU'; n = 1.5
else:
    marker = 'X'; color = '#bcbd22'; coordination = 'XX'

metal_rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}

nist = {
    'Ti': {'M': 1, 'O': 2, 'H_form': -944.747, 'G_form': -889.406},
    'V': {'M': 2, 'O': 5, 'H_form': -1550.590, 'G_form': -1419.359},
    'Cr': {'M': 2, 'O': 3, 'H_form': -1139.701, 'G_form': -1058.067},
    'Mn': {'M': 1, 'O': 1, 'H_form': -385.221, 'G_form': -362.898},
    'Fe': {'M': 2, 'O': 3, 'H_form': -824.248, 'G_form': -742.294},
    'Co': {'M': 3, 'O': 4, 'H_form': -910.020, 'G_form': -794.901},
    'Ni': {'M': 1, 'O': 1, 'H_form': -239.701, 'G_form': -211.539},
    'Cu': {'M': 1, 'O': 1, 'H_form': -156.063, 'G_form': -128.292},
}

root = '/pscratch/sd/j/jiuy97/3_V_bulk/metal/'
if os.path.exists(root):
    exp_path = '/pscratch/sd/j/jiuy97/3_V_bulk/oxide/monoxides.tsv'
    metal_path = '/pscratch/sd/j/jiuy97/3_V_bulk/metal/0_min/energy_norm_energy.tsv'
    oxide_path = '/pscratch/sd/j/jiuy97/3_V_bulk/oxide/0_min/energy_norm_energy.tsv'
    path = '/pscratch/sd/j/jiuy97/3_V_bulk/metal/merged_norm_energy.tsv'
else:
    exp_path = '/Users/jiuy97/Desktop/3_V_bulk/oxide/monoxides.tsv'
    metal_path = '/Users/jiuy97/Desktop/3_V_bulk/metal/0_min/energy_norm_energy.tsv'
    oxide_path = '/Users/jiuy97/Desktop/3_V_bulk/oxide/0_min/energy_norm_energy.tsv'
    path = '/Users/jiuy97/Desktop/3_V_bulk/metal/merged_norm_energy.tsv'

exp_df = pd.read_csv(exp_path, delimiter='\t')
metal_df = pd.read_csv(metal_path, delimiter='\t').iloc[:, 1:]
oxide_df = pd.read_csv(oxide_path, delimiter='\t').iloc[:, 1:]
df = pd.read_csv(path, delimiter='\t').iloc[:, 1:]

exp_df['dH_form'] = exp_df['dH_form'] / 96.48
metal_df.index = list(nist.keys())
oxide_df.index = list(nist.keys())
df.index = metal_rows['3d']

E_H2O = -14.23919983
E_H2 = -6.77409008

ZPE_H2O = 0.558
ZPE_H2 = 0.273

Cp_H2O = 0.10
Cp_H2 = 0.09

Ref_H2 = E_H2 + ZPE_H2 + Cp_H2
Ref_H2O = E_H2O + ZPE_H2O + Cp_H2O
Ref_O = Ref_H2O - Ref_H2 + 2.46 # 2.506
# print(Ref_H2O - Ref_H2)
# print(Ref_O)

for element, data in nist.items():
    data['H_form'] = data['H_form'] / data['M'] / 96.48
    data['OtoM'] = data['O'] / data['M']
    data['E_oxide'] = oxide_df.loc[element, 'energy']
    data['E_metal'] = data['E_oxide'] - data['H_form'] - data['OtoM'] * Ref_O

for metal in metal_rows['3d']:
    if metal in nist:
        df.loc[metal, '3d'] = nist[metal]['E_metal']

if os.path.exists(root):
    df.to_csv('/pscratch/sd/j/jiuy97/3_V_bulk/metal/corrected_norm_energy.tsv', sep='\t', float_format='%.4f')
else:
    df.to_csv('/Users/jiuy97/Desktop/3_V_bulk/metal/corrected_norm_energy.tsv', sep='\t', float_format='%.4f')

energy_path = './energy_norm_energy.tsv'
if not os.path.exists(energy_path):
    exit(1)

energy_df = pd.read_csv(energy_path, delimiter='\t', index_col=0)
formation = pd.DataFrame(index=energy_df.index, columns=energy_df.columns)

for row in metal_rows:
    if metal_rows[row] == energy_df.index.tolist():
        formation = energy_df.sub(df[row].values, axis=0) - n * Ref_O
        break

formation.to_csv(tsv_filename, sep='\t', float_format='%.4f')
print(f"Merged data saved to {tsv_filename}")

plt.figure(figsize=(8, 6), dpi=300)

for j, column in enumerate(formation.columns):
    x = formation.index
    y = formation[column]
    plt.plot(x, y, marker=marker, color=color, label='cal.')

for i, row in exp_df.iterrows():
    if row['row'] in metal_rows and row['Coordination'] == coordination:
        plt.scatter(row['numb'], row['dH_form'], marker=marker, color=color, 
                    edgecolors=color, facecolors='white', label='exp.')

plt.xlim(-0.5, len(x) - 0.5)
plt.xticks(np.arange(len(x)), x)
plt.xlabel('Metal (MO)')
plt.ylabel('Formation energy (eV/MO)')
plt.tight_layout()
plt.gcf().savefig(png_filename, bbox_inches="tight")
print(f"Figure saved as {png_filename}")
plt.close()
