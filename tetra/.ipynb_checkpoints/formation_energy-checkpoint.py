import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

png_filename = f"energy_norm_formation.png"   
tsv_filename = f"energy_norm_formation.tsv"

print(f"\033[92m{os.getcwd()}\033[0m")
if '1_Tetrahedral_WZ' in os.getcwd():
    marker = '>'; color = '#d62728'; coordination = 'WZ'
elif '2_Tetrahedral_ZB' in os.getcwd():
    marker = '<'; color = '#ff7f0e'; coordination = 'ZB'
elif '3_Tetragonal_LT' in os.getcwd():
    marker = 'o'; color = '#ffd70e'; coordination = 'LT'
elif '4_Square_Planar_TN' in os.getcwd():
    marker = 's'; color = '#2ca02c'; coordination = 'TN'
elif '5_Square_Planar_NB' in os.getcwd():
    marker = 'p'; color = '#279ff2'; coordination = 'NB'
elif '6_Octahedral_RS' in os.getcwd():
    marker = 'd'; color = '#9467bd'; coordination = 'RS'
else:
    marker = 'x'; color = 'k'; coordination = 'XX'

metal_rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
    }

nist = {
    'Ti': {'M': 1, 'O': 2, 'H_form': -944.747, 'G_form': -889.406}, # 1692 1809 mp-2657 Titanium Dioxide (Rutile)
    'V': {'M': 2, 'O': 5, 'H_form': -1550.590, 'G_form': -1419.359}, # 1780 1897 mp-25279 Divanadium Pentaoxide
    'Cr': {'M': 2, 'O': 3, 'H_form': -1139.701, 'G_form': -1058.067}, # 573 688 mp-19399 Dichromium Trioxide
    'Mn': {'M': 1, 'O': 1, 'H_form': -385.221, 'G_form': -362.898}, # 1046 1162 mp-19006 Manganese Oxide
    'Fe': {'M': 2, 'O': 3, 'H_form': -824.248, 'G_form': -742.294}, # 702 817 mp-19770 Hematite
    'Co': {'M': 3, 'O': 4, 'H_form': -910.020, 'G_form': -794.901}, # 544 659 mp-18748 Tricobalt Tetraoxide
    'Ni': {'M': 1, 'O': 1, 'H_form': -239.701, 'G_form': -211.539}, # 1213 1330 mp-19009 Nickel Oxide
    'Cu': {'M': 1, 'O': 1, 'H_form': -156.063, 'G_form': -128.292}, # 620 735 mp-704645 Copper Monoxide
    }

exp_path = '/pscratch/sd/j/jiuy97/3_V_bulk/oxide/monoxides.tsv'
metal_path = '/pscratch/sd/j/jiuy97/3_V_bulk/metal/0_min/energy_norm_energy.tsv'
oxide_path = '/pscratch/sd/j/jiuy97/3_V_bulk/oxide/0_min/energy_norm_energy.tsv'
path = '/pscratch/sd/j/jiuy97/3_V_bulk/metal/merged_norm_energy.tsv'

exp_df = pd.read_csv(exp_path, delimiter='\t')
metal_df = pd.read_csv(metal_path, delimiter='\t').iloc[:, 1:]
oxide_df = pd.read_csv(oxide_path, delimiter='\t').iloc[:, 1:]
df = pd.read_csv(path, delimiter='\t').iloc[:, 1:]

exp_df['dH_form'] = exp_df['dH_form'] / 96.48
metal_df.index = list(nist.keys())
oxide_df.index = list(nist.keys())
df.index = metal_rows['3d']

min_values = df.iloc[:, :3].min(axis=1)
df = df.iloc[:, 3:]

E_H2O = -14.23919983
E_H2 = -6.77409008

ZPE_H2O = 0.558
ZPE_H2 = 0.273

Cp_H2O = 0.10
Cp_H2 = 0.09

Ref_H2 = E_H2 + ZPE_H2 + Cp_H2
Ref_H2O = E_H2O + ZPE_H2O + Cp_H2O
Ref_O = Ref_H2O - Ref_H2 + 2.506

for element, data in nist.items():
    # data['G_form'] = data['G_form'] / data['M'] / 96.48
    # data['OtoM'] = data['O'] / data['M']
    # data['G_oxide'] = oxide_df.loc[element, 'energy'] - oxide_df.loc[element, 'TS'] + oxide_df.loc[element, 'ZPE']
    # data['G_metal'] = data['G_oxide'] - data['G_form'] - data['OtoM'] * G_oxygen
    # data['E_metal'] = data['G_metal'] + metal_df.loc[element, 'TS'] - metal_df.loc[element, 'ZPE']
    data['H_form'] = data['H_form'] / data['M'] / 96.48
    data['OtoM'] = data['O'] / data['M']
    data['E_oxide'] = oxide_df.loc[element, 'energy'] # - oxide_df.loc[element, 'TS'] + oxide_df.loc[element, 'ZPE']
    data['E_metal'] = data['E_oxide'] - data['H_form'] - data['OtoM'] * Ref_O

for i, metal in enumerate(metal_rows['3d']):
    if metal in nist:
        min_values.loc[metal] = nist[metal]['E_metal']
df.insert(0, '3d', min_values)
df.to_csv('/pscratch/sd/j/jiuy97/3_V_bulk/metal/corrected_norm_energy.tsv', sep='\t')

energy_path = './energy_norm_energy.tsv'
if not os.path.exists(energy_path):
    exit(1)
else:
    energy_df = pd.read_csv(energy_path, delimiter='\t', index_col=0)
    formation = pd.DataFrame(index=energy_df.index, columns=energy_df.columns)

for row in metal_rows:
    if metal_rows[row] == energy_df.index.tolist():
        formation = energy_df.sub(df[row].values, axis=0) - Ref_O # G_oxygen
        break
formation.to_csv(tsv_filename, sep='\t')
print(f"Merged data saved to {tsv_filename}")

plt.figure(figsize=(8, 6))
# for j, column in enumerate(formation.columns):
#     filtered_x = []
#     filtered_values = []
#     x = formation.index
#     values = formation[column]
#     for i, v in enumerate(values):
#         if not np.isnan(v):
#             filtered_x.append(i)
#             filtered_values.append(v)
#     if not filtered_values:
#         print(f"No values found for pattern: {column}")
#         continue
#     plt.plot(filtered_x, filtered_values, marker=marker, color=color, label='cal.')
    
for j, column in enumerate(formation.columns):
    x = formation.index
    y = formation[column]
    plt.plot(x, y, marker=marker, color=color, label='cal.')
    
for i in exp_df.index:
    if exp_df['row'][i] == row and exp_df['Coordination'][i] == coordination:
        plt.scatter(exp_df['numb'][i], exp_df['dH_form'][i], 
                    marker=marker, color=color, edgecolors=color, facecolors='white', label='exp.')

plt.xlim(-0.5, len(x)-0.5)
plt.xticks(np.arange(len(x)), x)
plt.xlabel('Metal (MO)')
plt.ylabel('Formation energy (eV/MO)')
plt.tight_layout()
plt.gcf().savefig(png_filename, bbox_inches="tight")
print(f"Figure saved as {png_filename}")
plt.close()