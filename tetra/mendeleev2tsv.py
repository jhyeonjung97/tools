import matplotlib.pyplot as plt
from mendeleev import element
import pandas as pd
import numpy as np

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
    for metal in metals[row]:
        elem = element(metal)
        for pattern in patterns:
            name = patterns[pattern]
            if 'ionenergies' in pattern:
                ion_index = int(pattern.split('[')[1].strip(']'))
                df.loc[metal, name] = elem.ionenergies.get(ion_index, np.nan)
            else:
                value = getattr(elem, pattern)
                if isinstance(value, float):
                    df.at[metal, name] = value
                else:
                    print(metal, name, value)
                    df.at[metal, name] = np.nan
                
print(df)
            
    # # if pattern == 'boiling_point' or pattern == 'melting_point':
    # if pattern == 'melting_point':
    #     for i in range(n):
    #         j = 13 * i + 12
    #         df.loc[j, '4d'] = 286.35
            
    # if pattern == 'evaporation_heat':
    #     for i in range(n):
    #         j = 13 * i + 6
    #         df.loc[j, '4d'] = 619.00

    # df.index = index_pattern

    # if n == 1:
    #     tsv_filename=f'mendeleev_{filename}.tsv'
    #     png_filename=f'mendeleev_{filename}.png'
    # else:
    #     tsv_filename=f'concat_{filename}.tsv'
    #     png_filename=f'concat_{filename}.png'
    
    # df.to_csv(f'{tsv_filename}', sep='\t', index=True, float_format='%.4f')

    # print(f"DataFrame saved as {tsv_filename}")

    # if n == 1:
    #     plt.figure()
    #     plt.plot(df.index, df['3d'], marker='d', color=colors[0], label='3d')
    #     plt.plot(df.index, df['4d'], marker='d', color=colors[1], label='4d')
    #     plt.plot(df.index, df['5d'], marker='d', color=colors[2], label='5d')
    #     plt.xticks(np.arange(len(indice)), indice)
    #     plt.xlabel('Metal (MO)')
    #     plt.ylabel(pattern.replace('_', ' ').title())
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(f'{png_filename}', bbox_inches="tight")
    #     print(f"Figure saved as {png_filename}")
    #     plt.close()