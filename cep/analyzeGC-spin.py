import matplotlib.pyplot as plt
import numpy as np
import os

# Folder names, labels, and colors for each spin state
folders = ['hs', 'is', 'ls']
labels = ['HS', 'IS', 'LS']
colors = ['tab:purple', 'tab:blue', 'tab:orange']

plt.figure(figsize=(8, 5))

# Read data from each folder and plot with quadratic fitting
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'GCFE_data_FULL.dat')
    data = np.loadtxt(filepath)
    potentials = data[:, 0]
    energies = data[:, 1]

    # Scatter plot of data points
    plt.scatter(potentials, energies, color=color, label=label)

    # Perform quadratic fitting
    coeffs = np.polyfit(potentials, energies, deg=2)
    poly = np.poly1d(coeffs)

    # Generate x values for smooth fitting curve
    x_fit = np.linspace(min(potentials), max(potentials), 200)
    y_fit = poly(x_fit)

    # Plot the fitted quadratic curve
    plt.plot(x_fit, y_fit, color=color, linestyle='--')

plt.xlabel("Potential vs. SHE")
plt.ylabel("Free Energy, eV")
plt.legend()
plt.tight_layout()
plt.savefig('spin.png')
plt.show()
plt.close()