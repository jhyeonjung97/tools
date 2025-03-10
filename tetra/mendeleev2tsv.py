import matplotlib.pyplot as plt
from mendeleev import element
import pandas as pd
import numpy as np

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}
patterns = [
    'group_id', 'atomic_number', 'atomic_volume', 'boiling_point', 'melting_point',
    'mass', 'density', 'dipole_polarizability', 'en_pauling', 'covalent_radius',
    'metallic_radius', 'vdw_radius', 'evaporation_heat', 'fusion_heat',
    'heat_of_formation', 'ionenergies[1]', 'ionenergies[2]', 'ionenergies[3]'
]

def get_data(element_symbol, atomic_property):
    try:
        elem = element(element_symbol)
        if 'ionenergies' in atomic_property:
            ion_index = int(atomic_property.split('[')[1].strip(']'))
            return elem.ionenergies.get(ion_index, None)  # Prevent KeyError
        else:
            return getattr(elem, atomic_property, None)
    except AttributeError:
        return None
    except (KeyError, IndexError):
        return None

for pattern in patterns:
    data = {
        '3d': [get_data(e, pattern) for e in metals['3d']],
        '4d': [get_data(e, pattern) for e in metals['4d']],
        '5d': [get_data(e, pattern) for e in metals['5d']],
    }

    df = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce')  # Convert all to numeric safely
    print(f"\nProperty: {pattern}")
    print(df.head())  # Show first few rows for debugging

            
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