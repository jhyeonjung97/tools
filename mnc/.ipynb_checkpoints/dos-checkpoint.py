import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def gaussian_smearing(x, y, sigma, num_points=1000):
    """Apply Gaussian smearing to DOS."""
    x_grid = np.linspace(x.min(), x.max(), num_points)
    smeared_y = np.zeros_like(x_grid)
    for xi, yi in zip(x, y):
        smeared_y += yi * np.exp(-0.5 * ((x_grid - xi) / sigma) ** 2)
    smeared_y /= (sigma * np.sqrt(2 * np.pi))  # Normalize
    return x_grid, smeared_y

def get_fermi_level(doscar_path="DOSCAR"):
    """Extract Fermi level from the DOSCAR file."""
    with open(doscar_path, 'r') as file:
        lines = file.readlines()
    fermi_level = float(lines[5].split()[3])  # Fermi level is the 4th value on the 6th line
    return fermi_level

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Plot projected DOS with Gaussian smearing and Fermi level alignment.")
parser.add_argument("--file", type=str, required=True, help="Path to the input DOS data file.")
parser.add_argument("--gaussian", type=float, default=0.05, help="Gaussian smearing value in eV (default: 0.05).")
parser.add_argument("--xrange", type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=(-8, 6),
                    help="Energy range for the x-axis (default: -8 to 6).")
args = parser.parse_args()

# Extract Fermi level from DOSCAR
fermi_level = get_fermi_level()

# Define column names
columns = [
    "energy", "s(up)", "s(down)", "p(up)", "p(down)",
    "dxy(up)", "dxy(down)", "dyz(up)", "dyz(down)",
    "dz2(up)", "dz2(down)", "dxz(up)", "dxz(down)",
    "dx2(up)", "dx2(down)"
]

# Read the DOS data
data = pd.read_csv(args.file, delim_whitespace=True, comment='#', names=columns)

# Align energy with the Fermi level
energy = data["energy"] - fermi_level

# Define orbitals to plot
orbitals = ["dxy", "dyz", "dz2", "dxz", "dx2"]

# Create subplots
fig, axes = plt.subplots(len(orbitals), 1, figsize=(8, len(orbitals) * 3), sharex=True)

for i, orbital in enumerate(orbitals):
    orbital_up = data[f"{orbital}(up)"]
    orbital_down = data[f"{orbital}(down)"]

    # Apply Gaussian smearing if sigma > 0
    if args.gaussian > 0:
        energy, orbital_up = gaussian_smearing(energy, orbital_up, args.gaussian)
        _, orbital_down = gaussian_smearing(energy, orbital_down, args.gaussian)

    # Plot on the respective subplot
    axes[i].plot(energy, orbital_up, label=f"{orbital} (up)", linestyle='-', linewidth=1.5)
    axes[i].plot(energy, orbital_down, label=f"{orbital} (down)", linestyle='--', linewidth=1.5)
    axes[i].axhline(0, color='black', linewidth=0.8, linestyle='--')
    axes[i].legend()
    axes[i].set_ylabel("DOS")
    axes[i].grid(alpha=0.3)
    axes[i].set_title(f"Projected DOS for {orbital} Orbital")

# Set shared x-axis labels
axes[-1].set_xlabel("Energy (eV) (relative to Fermi level)")
plt.xlim(args.xrange)

# Adjust layout
plt.tight_layout()
plt.show()
