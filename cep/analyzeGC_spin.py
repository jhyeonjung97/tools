import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Folder names, labels, and colors for each spin state
folders = ['hs', 'is', 'ls']
labels = ['HS', 'IS', 'LS']
colors = ['tab:orange', 'tab:blue', 'tab:green']

# Check which folders exist
existing_folders = [f for f in folders if os.path.exists(f)]
if not existing_folders:
    print("Error: No spin state folders found.")
    exit(1)

# Update folders, labels, and colors to only include existing ones (use original list for index)
labels = [labels[folders.index(f)] for f in existing_folders]
colors = [colors[folders.index(f)] for f in existing_folders]
folders = existing_folders

# Plot GC Free Energy
plt.figure(figsize=(4, 3))

# Loop through each folder and plot if the data file exists
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'gc_energy.dat')
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
        plt.scatter(potentials, energies, color=color, zorder=10)

        # Perform quadratic fitting
        coeffs = np.polyfit(potentials, energies, deg=2)
        poly = np.poly1d(coeffs)
        
        # Calculate R²
        y_pred = poly(potentials)
        ss_res = np.sum((energies - y_pred) ** 2)
        ss_tot = np.sum((energies - np.mean(energies)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Format equation: y = ax² + bx + c
        a, b, c = coeffs[0], coeffs[1], coeffs[2]
        b_sign = '' if b < 0 else '+'
        c_sign = '' if c < 0 else '+'
        eq_text = r'{}: $y = {:.2f}x^2 {}{:.2f}x {}{:.2f}$ (R² = {:.2f})'.format(
            label, a, b_sign, b, c_sign, c, r2
        )

        # Generate x values for smooth fitting curve
        x_fit = np.linspace(min(potentials), max(potentials), 200)
        y_fit = poly(x_fit)

        # Plot the fitted quadratic curve
        plt.plot(x_fit, y_fit, color=color)
        plt.plot([], [], marker='o', color=color, label=eq_text)
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        continue

plt.xlabel("Potential vs. SHE", fontsize=12)
plt.ylabel("Free Energy, eV", fontsize=12)
plt.legend(fontsize=10, loc='lower right', bbox_to_anchor=(1.03, 1.0), ncol=1, labelspacing=0.2)
plt.savefig('spin_full.png', dpi=300, bbox_inches='tight')
plt.xlim(0.0, 1.0)
plt.savefig('spin.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot dz values
plt.figure(figsize=(4, 3))

# Loop through each folder and plot if the data file exists
for folder, label, color in zip(folders, labels, colors):
    filepath = os.path.join(folder, 'dz.dat')
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
        
        # Calculate R²
        y_pred = poly(potentials)
        ss_res = np.sum((dz_values - y_pred) ** 2)
        ss_tot = np.sum((dz_values - np.mean(dz_values)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Format equation: y = ax² + bx + c
        a, b, c = coeffs[0], coeffs[1], coeffs[2]
        b_sign = '' if b < 0 else '+'
        c_sign = '' if c < 0 else '+'
        eq_text = r'{}: $y = {:.2f}x^2 {}{:.2f}x {}{:.2f}$ (R² = {:.2f})'.format(
            label, a, b_sign, abs(b), c_sign, abs(c), r2
        )

        # Generate x values for smooth fitting curve
        x_fit = np.linspace(min(potentials), max(potentials), 200)
        y_fit = poly(x_fit)

        # Plot the fitted quadratic curve
        plt.plot(x_fit, y_fit, color=color)
        plt.plot([], [], marker='o', color=color, label=eq_text)
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        continue

plt.xlabel("Potential vs. SHE", fontsize=12)
plt.ylabel("dz (Å)", fontsize=12)
plt.ylim(0, 1.0)
plt.legend(fontsize=10, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=1, labelspacing=0.5)
plt.savefig('spin_dz_full.png', dpi=300, bbox_inches='tight')
plt.xlim(0.0, 1.0)
plt.savefig('spin_dz.png', dpi=300, bbox_inches='tight')
plt.close()
