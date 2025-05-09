import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dir_now = os.getcwd()
print(f"\033[92m{dir_now}\033[0m")

# Read input data
Ef_oxide = pd.read_csv('./energy_norm_formation.tsv', delimiter='\t', index_col=0)
sublimation_heat_path = '/pscratch/sd/j/jiuy97/3_V_bulk/6_Octahedral_RS/mendeleev_sublimation_heat.tsv'
if not os.path.exists(sublimation_heat_path):
    sublimation_heat_path = '/Users/jiuy97/Desktop/3_V_bulk/6_Octahedral_RS/mendeleev_sublimation_heat.tsv'
if not os.path.exists(sublimation_heat_path):
    sublimation_heat_path = '/Users/hailey/Desktop/3_V_bulk/6_Octahedral_RS/mendeleev_sublimation_heat.tsv'
Ec_metal = pd.read_csv(sublimation_heat_path, delimiter='\t', index_col=0)

# Define the rows for different series of metals
rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}

# Fixed value for oxygen cohesive energy
Ec_oxygen = 5.1614  # eV

# Initialize an empty DataFrame for cohesive energy calculations
Ec_oxide = pd.DataFrame(index=Ef_oxide.index)
Ec_metal.index = Ef_oxide.index

# Determine the current working directory and perform calculations accordingly

if 'fm' in dir_now or '3d' in dir_now:
    row = '3d'
elif '4d' in dir_now:
    row = '4d'
elif '5d' in dir_now:
    row = '5d'
else:
    raise ValueError("The current directory does not match any expected pattern ('1_afm', '2_fm', '3d', '4d', '5d').")

# Calculate the cohesive energy
Ec_oxide['energy'] = Ec_metal[row] / 96.48 + Ec_oxygen / 2 - Ef_oxide['energy']

# Save the calculated cohesive energies to a TSV file
Ec_oxide.to_csv('energy_norm_cohesive.tsv', sep='\t', index=True, float_format='%.4f')
print(f"Data saved to energy_norm_cohesive.tsv")

# Set plotting parameters based on the current working directory
if '1_Tetrahedral_WZ' in dir_now:
    marker = '>'; color = '#d62728'
elif '2_Tetrahedral_ZB' in dir_now:
    marker = '<'; color = '#ff7f0e'
elif '3_SquarePlanar_TN' in dir_now:
    marker = 'o'; color = '#ffd70e'
elif '4_SquarePlanar_PD' in dir_now:
    marker = 's'; color = '#2ca02c'
elif '5_SquarePlanar_NB' in dir_now:
    marker = 'p'; color = '#17becf'
elif '6_Octahedral_RS' in dir_now:
    marker = 'd'; color = '#9467bd'
elif '7_Pyramidal_LT' in dir_now:
    marker = 'h'; color = '#8c564b'
elif '8_Tetrahedral_AQ' in dir_now:
    marker = '^'; color = '#e377c2'
elif '9_SquarePlanar_AU' in dir_now:
    marker = 'v'; color = '#7f7f7f'
else:
    marker = 'X'; color = '#bcbd22'

# Plotting the cohesive energy
plt.figure(figsize=(8, 6), dpi=300)
# x = []
# filtered_values = []
# for i, v in enumerate(Ec_oxide['energy']):
#     if not np.isnan(v): 
#         x.append(i)
#         filtered_values.append(v)
# plt.plot(x, filtered_values, marker=marker, color=color)
x = range(len(Ec_oxide['energy']))
y = Ec_oxide['energy']
plt.plot(x, y, marker=marker, color=color)
plt.xlabel('Metal (MO)')
plt.ylabel('Cohesive energy (eV)')
plt.xticks(np.arange(len(Ec_oxide.index)), Ec_oxide.index)
plt.xlim(-0.5, len(Ec_oxide.index))  # Add space before the first x-tick
plt.tight_layout()
plt.savefig('energy_norm_cohesive.png', bbox_inches="tight")
plt.close()
print("Figure saved as energy_norm_cohesive.png")