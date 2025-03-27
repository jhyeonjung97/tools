import os
import numpy as np
import pandas as pd
from mendeleev import element
import matplotlib.pyplot as plt

root = '/Users/hailey/Desktop/7_V_bulk'

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])]
colors = plt.cm.Purples(np.linspace(0.4, 0.9, 3))
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

# 데이터프레임 초기화 시 데이터 타입 지정
df = pd.DataFrame(index=sum(metals.values(), []), dtype='object')
df['row'] = pd.Series(dtype='object')
df['numb'] = pd.Series(dtype='int')

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
            elif metal == 'Sn' and pattern in ['boiling_point', 'melting_point']:
                # Sn의 경우 특별 처리
                if pattern == 'boiling_point':
                    df.at[metal, column] = 2875  # Sn의 끓는점 (K)
                elif pattern == 'melting_point':
                    df.at[metal, column] = 505.08  # Sn의 녹는점 (K)
            else:
                value = getattr(elem, pattern)
                if isinstance(value, dict):
                    df.at[metal, column] = value.get('gray', value)
                else:
                    df.at[metal, column] = value

df['ion12'] = df['ion1'] + df['ion2']
save_path = os.path.join(root, 'figures')
df.to_csv(f'{save_path}/mendeleev_data.csv', sep=',')
df.to_csv(f'{save_path}/mendeleev_data.tsv', sep='\t', float_format='%.2f')
print(f"Data saved as mendeleev_data.csv")

for column in df.columns:
    if column in ['row', 'numb']:
        continue
    elif column.startswith('ion'):
        pattern = f'ionenergies{column.replace("ion", "")}'
    else:
        pattern = next((key for key, value in patterns.items() if value == column), None)
    pngname = f'mendeleev_{pattern}.png'
        
    plt.figure(figsize=(8, 6), dpi=100)
    plt.plot(df[df['row'] == '3d']['numb'], df[df['row'] == '3d'][column], 
             marker='d', color=colors[0], label='3d')
    plt.plot(df[df['row'] == '4d']['numb'], df[df['row'] == '4d'][column], 
             marker='d', color=colors[1], label='4d')
    plt.plot(df[df['row'] == '5d']['numb'], df[df['row'] == '5d'][column], 
             marker='d', color=colors[2], label='5d')
    plt.xticks(np.arange(len(indice)), indice)
    plt.xlabel('Metal (MO)')
    plt.ylabel(pattern.replace('_', ' ').title())
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{save_path}/{pngname}', bbox_inches="tight")
    print(f"Figure saved as {pngname}")
    plt.close()