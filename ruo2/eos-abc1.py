import os
import pandas as pd
from ase.io import read
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

base_path = "~/Desktop/3_RuO2/4_high_valence/1_M-RuO2/4_Re"
base_path = os.path.expanduser(base_path)
dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
dirs.sort()

figsize = (6.4, 4.8)
plt.figure(figsize=figsize)
df = pd.DataFrame(columns=['a', 'b', 'c', 'energy'])
for dir in dirs:
    parts = dir.split('_')
    a = float(parts[1])
    b = float(parts[2])
    c = float(parts[3])
    json_path = os.path.join(base_path, dir, 'final_with_calculator.json')
    if not os.path.exists(json_path):
        continue
    atoms = read(json_path)
    energy = atoms.get_potential_energy()
    df.loc[len(df)] = [a, b, c, energy]
vmin = df['energy'].min()
df['energy'] = (df['energy'] - vmin)
df.to_csv('eos-abc.csv', index=False)
df = df.loc[df.groupby(['a', 'c'])['energy'].idxmin()].reset_index(drop=True)
print(df)
plt.scatter(df['a'], df['c'], c=df['energy'], edgecolors='black', linewidths=0.5, cmap='RdBu', zorder=10)
plt.tricontourf(df['a'], df['c'], df['energy'], levels=30, cmap='RdBu', zorder=0)
plt.colorbar(label='Relative energy (eV)')
plt.gca().set_aspect('equal', adjustable='box')
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.01))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.xlabel('Lattice constant a (Å)')
plt.ylabel('Lattice constant c (Å)')
plt.tight_layout()
plt.savefig('eos-ac.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()