import matplotlib.pyplot as plt
from mendeleev import element
import pandas as pd
import numpy as np
import argparse
import os

print(f"\033[92m{os.getcwd()}\033[0m")
parser = argparse.ArgumentParser(description='Generate TSV files with specified atomic properties and plot if n=1.')
parser.add_argument('-p', '--patterns', required=True, nargs='+', help='List of atomic properties to retrieve.')
parser.add_argument('-n', '--number', default=1, type=int, help='Number of repeat')

args = parser.parse_args()
n = args.number
m = n * 13

metal_rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
    }
indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metal_rows['3d'], metal_rows['4d'], metal_rows['5d'])]
colors = plt.cm.Purples(np.linspace(0.4, 0.9, 3))

def get_data(element_symbol, atomic_property):
    try:
        elem = element(element_symbol)
        if 'ionenergies' in atomic_property:
            ion_index = int(atomic_property.split('[')[1].strip(']'))
            return elem.ionenergies[ion_index]
        else:
            return getattr(elem, atomic_property)
    except AttributeError:
        return None
    except (KeyError, IndexError):
        return None

index_pattern = list(range(13)) * n
index_pattern = index_pattern[:m]

for pattern in args.patterns:
    if 'ionenergies' in pattern:
        ion_index = int(pattern.split('[')[1].strip(']'))
        filename = f'ionenergies_{ion_index}'
    else:
        filename = pattern

    data = {
        '3d': [get_data(e, pattern) for e in metal_rows['3d']] * n,
        '4d': [get_data(e, pattern) for e in metal_rows['4d']] * n,
        '5d': [get_data(e, pattern) for e in metal_rows['5d']] * n
    }

    df = pd.DataFrame(data).iloc[:m]
    
    if pattern == 'boiling_point' or pattern == 'melting_point':
        for i in range(n):
            j = 13 * i + 12
            df.loc[j, '4d'] = df.loc[j, '4d']['gray']
    
    if pattern == 'evaporation_heat':
        for i in range(n):
            j = 13 * i + 6
            df.loc[j, '4d'] = 619.0

    df.index = index_pattern

    if n == 1:
        tsv_filename=f'mendeleev_{filename}.tsv'
        png_filename=f'mendeleev_{filename}.png'
    else:
        tsv_filename=f'concat_{filename}.tsv'
        png_filename=f'concat_{filename}.png'
    
    df.to_csv(f'{tsv_filename}', sep='\t', index=True, float_format='%.2f')

    print(f"DataFrame saved as {tsv_filename}")

    if n == 1:
        plt.figure()
        plt.plot(df.index, df['3d'], marker='d', color=colors[0], label='3d')
        plt.plot(df.index, df['4d'], marker='d', color=colors[1], label='4d')
        plt.plot(df.index, df['5d'], marker='d', color=colors[2], label='5d')
        plt.xticks(np.arange(len(indice)), indice)
        plt.xlabel('Metal (MO)')
        plt.ylabel(pattern.replace('_', ' ').title())
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{png_filename}', bbox_inches="tight")
        print(f"Figure saved as {png_filename}")
        plt.close()