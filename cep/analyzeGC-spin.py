import matplotlib.pyplot as plt
import numpy as np
import os

# 폴더 이름과 레이블 지정
folders = ['hs', 'is', 'ls']
labels = ['High Spin', 'Intermediate Spin', 'Low Spin']
colors = ['red', 'green', 'blue']

plt.figure(figsize=(8, 6))

# 각 폴더에서 데이터 읽고 플롯
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'GCFE_data_FULL.dat')
    data = np.loadtxt(filepath)
    potentials = data[:, 0]
    energies = data[:, 1]
    plt.plot(potentials, energies, label=label, color=color)

plt.xlabel('Applied Potential (V)')
plt.ylabel('Gibbs Free Energy (eV)')
plt.title('Gibbs Free Energy vs. Applied Potential')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
