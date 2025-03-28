import matplotlib.pyplot as plt
import numpy as np
import os

# Folder names, labels, and colors for each spin state
folders = ['hs', 'is', 'ls']
labels = ['HS', 'IS', 'LS']
colors = ['tab:orange', 'tab:blue', 'tab:green']

# Check which folders exist
existing_folders = [f for f in folders if os.path.exists(f)]
if not existing_folders:
    print("Error: No spin state folders found.")
    exit(1)

# Update folders, labels, and colors to only include existing ones
folders = existing_folders
labels = [labels[folders.index(f)] for f in existing_folders]
colors = [colors[folders.index(f)] for f in existing_folders]

# Plot GC Free Energy
plt.figure(figsize=(4, 3))

# Loop through each folder and plot if the data file exists
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'GCFE_data_FULL.dat')
    if not os.path.exists(filepath):
        print(f"Warning: File not found: {filepath}, skipping.")
        continue

    try:
        data = np.loadtxt(filepath)
        if len(data) == 0:
            print(f"Warning: Empty file: {filepath}, skipping.")
            continue
            
        potentials = data[:, 0]
        energies = data[:, 1]

        # Scatter plot of data points
        plt.scatter(potentials, energies, color=color)

        # Perform quadratic fitting
        coeffs = np.polyfit(potentials, energies, deg=2)
        poly = np.poly1d(coeffs)

        # Generate x values for smooth fitting curve
        x_fit = np.linspace(min(potentials), max(potentials), 200)
        y_fit = poly(x_fit)

        # Plot the fitted quadratic curve
        plt.plot(x_fit, y_fit, color=color)
        plt.plot([], [], marker='o', color=color, label=label)
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        continue

plt.xlabel("Potential vs. SHE")
plt.ylabel("Free Energy, eV")
plt.legend()
plt.tight_layout()
plt.savefig('spin_full.png')
plt.xlim(0.0, 1.0)
plt.savefig('spin.png')
plt.close()

# Plot dz values
plt.figure(figsize=(4, 3))

# Loop through each folder and plot if the data file exists
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'dz_data_FULL.dat')
    if not os.path.exists(filepath):
        print(f"Warning: File not found: {filepath}, skipping.")
        continue

    try:
        data = np.loadtxt(filepath)
        if len(data) == 0:
            print(f"Warning: Empty file: {filepath}, skipping.")
            continue
            
        potentials = data[:, 0]
        dz_values = data[:, 1]

        # Scatter plot of data points
        plt.scatter(potentials, dz_values, color=color)

        # Perform quadratic fitting
        coeffs = np.polyfit(potentials, dz_values, deg=2)
        poly = np.poly1d(coeffs)

        # Generate x values for smooth fitting curve
        x_fit = np.linspace(min(potentials), max(potentials), 200)
        y_fit = poly(x_fit)

        # Plot the fitted quadratic curve
        plt.plot(x_fit, y_fit, color=color)
        plt.plot([], [], marker='o', color=color, label=label)
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        continue

plt.xlabel("Potential vs. SHE")
plt.ylabel("dz (Å)")
plt.ylim(0, 1.0)
plt.legend()
plt.tight_layout()
plt.savefig('spin_dz_full.png')
plt.xlim(0.0, 1.0)
plt.savefig('spin_dz.png')
plt.close()
