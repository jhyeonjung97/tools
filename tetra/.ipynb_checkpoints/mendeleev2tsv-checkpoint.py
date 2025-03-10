import os
import numpy as np
import pandas as pd
from mendeleev import element
import matplotlib.pyplot as plt

root = '/pscratch/sd/j/jiuy97/7_V_bulk'

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}
patterns = {
    'group_id': 'group', 
    'atomic_number': 'Natom', 
    'atomic_volume': 'Vatom', 
    'boiling_point': 'Tboil', 
    'melting_point': 'Tmelt',
    'mass': 'mass', 
    'density': 'density', 
    'dipole_polarizability': 'dipole', 
    'en_pauling': 'pauling', 
    'covalent_radius': 'Rcoval',
    'metallic_radius': 'Rmetal', 
    'vdw_radius': 'Rvdw', 
    'evaporation_heat': 'Hevap', 
    'fusion_heat': 'Hfus',
    'heat_of_formation': 'Hform', 
    'ionenergies[1]': 'ion1', 
    'ionenergies[2]': 'ion2', 
    'ionenergies[3]': 'ion3',
}

df = pd.DataFrame()
for row in metals.keys():
    for m, metal in enumerate(metals[row]):
        df.at[metal, 'row'] = row
        df.at[metal, 'numb'] = m
        elem = element(metal)
        for pattern in patterns:
            column = patterns[pattern]
            if 'ionenergies' in pattern:
                ion_index = int(pattern.split('[')[1].strip(']'))
                df.at[metal, column] = elem.ionenergies.get(ion_index, np.nan)
            elif metal == 'Sn' and (pattern == 'boiling_point' or pattern == 'melting_point'):
                df.at[metal, column] = getattr(elem, pattern)['gray']
            else:
                df.at[metal, column] = getattr(elem, pattern)

df['ion12'] = df['ion1'] + df['ion2']
save_path = os.path.join(root, 'figures')
df.to_csv(f'{save_path}/mendeleev_data.csv', sep=',')
df.to_csv(f'{save_path}/mendeleev_data.tsv', sep='\t')

for column in df.columns:
    pngname = f'mendeleev_{column}.png'
        
    plt.figure()
    plt.plot(df[df['row'] == '3d']['numb'], df[df['row'] == '3d'][column], label='3d')
    plt.plot(df[df['row'] == '4d']['numb'], df[df['row'] == '4d'][column], label='4d')
    plt.plot(df[df['row'] == '5d']['numb'], df[df['row'] == '5d'][column], label='5d')
    plt.xlabel('Metal (MO)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{png_filename}', bbox_inches="tight")
    print(f"Figure saved as {png_filename}")
    plt.close()