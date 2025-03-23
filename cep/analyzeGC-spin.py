import matplotlib.pyplot as plt
import numpy as np
import os

# 폴더 이름과 레이블 지정
folders = ['hs', 'is', 'ls']
labels = ['HS', 'IS', 'LS']
colors = ['purple', 'blue', 'orange']

plt.figure(figsize=(8, 5))

# 각 폴더에서 데이터 읽고 플롯
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'GCFE_data_FULL.dat')
    data = np.loadtxt(filepath)
    potentials = data[:, 0]
    energies = data[:, 1]
    plt.plot(potentials, energies, label=label, color=color)

plt.xlabel("Potential vs. SHE")
plt.ylabel("Free Energy, eV")
plt.legend()
plt.tight_layout()
plt.savefig('spin.png')
plt.show()
plt.close()