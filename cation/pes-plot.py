import os
import numpy as np
import pandas as pd
from ase.io import read
import matplotlib.pyplot as plt

root = '/Users/hailey/Desktop/1_cation/3_PES'

df_energy = pd.DataFrame(index=range(31))
df_charge = pd.DataFrame(index=range(31))
for solv in ['1_cation_water', '2_cation_hydroxide']:
    solv_name = 'water' if solv == '1_cation_water' else 'hydro'
    for dir in ['pes', 'relaxed-pes']:
        dir_name = 'fixed' if dir == 'pes' else 'relaxed'
        for subdir in ['0_H', '1_Li', '2_Na', '3_K', '4_Rb', '5_Cs']:
            subdir_name = subdir.split('_')[1]
            for i in range(31):
                file_energy = os.path.join(root, solv, dir, subdir, f'{i:02d}_', 'final_with_calculator.json')
                file_charge = os.path.join(root, solv, dir, subdir, f'{i:02d}_', 'atoms_with_calculator.json')
                if os.path.isfile(file_energy):
                    energy = read(file_energy).get_potential_energy()
                    df_energy.loc[i, f'{solv_name}_{dir_name}_{subdir_name}'] = energy
                    # print(solv_name, dir_name, subdir_name, energy)

df_energy = df_energy.sub(df_energy.loc[0], axis=1)
print(df_energy.filter(like='water_relaxed').loc[0:4])
print(df_energy.filter(like='hydro_relaxed').loc[0:4])

plt.figure(figsize=(10, 6))
for col in df_energy.columns:
    if 'H' in col:
        color = 'black'
    elif 'Li' in col:
        color = 'tab:red'
    elif 'Na' in col:
        color = 'tab:orange'
    elif 'K' in col:
        color = 'tab:green'
    elif 'Rb' in col:
        color = 'tab:blue'
    elif 'Cs' in col:
        color = 'tab:purple'
    edgecolor = color

    if 'fixed' in col:
        facecolor = 'white'
    elif 'relaxed' in col:
        facecolor = edgecolor

    if 'water' in col:
        marker = 'o'
    elif 'hydro' in col:
        marker = 's'

    plt.plot(df_energy.index, df_energy[col], marker=marker, label=col, linestyle='-',
             color=color, markerfacecolor=facecolor, markeredgecolor=edgecolor, markersize=5)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(root, 'pes-energy.png'), bbox_inches="tight")
plt.tight_layout()
plt.show()