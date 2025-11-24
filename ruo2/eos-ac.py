import os
import pandas as pd
from ase.io import read
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

base_path = "~/Desktop/3_RuO2/4_high_valence/1_M-RuO2/0_Ru"
base_path = os.path.expanduser(base_path)
dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
dirs.sort()

figsize = (6.4/6.4*5.4, 4.8/6.4*5.4)
plt.figure(figsize=figsize)
df = pd.DataFrame(columns=['a', 'b', 'energy'])
for dir in dirs:
    parts = dir.split('_')
    a = float(parts[1])
    b = float(parts[2])
    atoms = read(os.path.join(base_path, dir, 'final_with_calculator.json'))
    energy = atoms.get_potential_energy()
    df.loc[len(df)] = [a, b, energy]
vmin = df['energy'].min()
df['energy'] = (df['energy'] - vmin)
df.to_csv('eos-ac.csv', index=False)
plt.scatter(df['a'], df['b'], c=df['energy'], edgecolors='black', linewidths=0.5, cmap='RdBu', zorder=10)
plt.tricontourf(df['a'], df['b'], df['energy'], levels=30, cmap='RdBu', zorder=0)
plt.colorbar(label='Relative energy (eV)')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Lattice constant a=b (Å)')
plt.ylabel('Lattice constant c (Å)')

# Tick 설정: 간격 0.01, 소수점 두 자리
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.01))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x:.2f}'))
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x:.2f}'))
plt.tight_layout()
plt.savefig('eos-ac.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()