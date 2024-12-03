import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def gaussian_smearing(x, y, sigma, num_points=1000):
    """Apply Gaussian smearing to DOS."""
    # Create a fine grid for the smeared DOS
    x_grid = np.linspace(x.min(), x.max(), num_points)
    smeared_y = np.zeros_like(x_grid)
    for xi, yi in zip(x, y):
        smeared_y += yi * np.exp(-0.5 * ((x_grid - xi) / sigma) ** 2)
    smeared_y /= (sigma * np.sqrt(2 * np.pi))  # Normalize
    return x_grid, smeared_y

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Plot projected DOS with Gaussian smearing.")
parser.add_argument("--file", type=str, required=True, help="Path to the input DOS data file.")
parser.add_argument("--gaussian", type=float, default=0.0, help="Gaussian smearing value in eV (default: 0, no smearing).")
args = parser.parse_args()

# Define column names
columns = [
    "energy", "s(up)", "s(down)", "p(up)", "p(down)",
    "dxy(up)", "dxy(down)", "dyz(up)", "dyz(down)",
    "dz2(up)", "dz2(down)", "dxz(up)", "dxz(down)",
    "dx2(up)", "dx2(down)"
]

# Read the data
data = pd.read_csv(args.file, delim_whitespace=True, comment='#', names=columns)

# Extract energy and dxy orbital data
energy = data["energy"]
dxy_up = data["dxy(up)"]
dxy_down = data["dxy(down)"]

# Apply Gaussian smearing if sigma > 0
if args.gaussian > 0:
    energy, dxy_up = gaussian_smearing(energy, dxy_up, args.gaussian)
    _, dxy_down = gaussian_smearing(data["energy"], dxy_down, args.gaussian)

# Plot the dxy orbital
plt.figure(figsize=(8, 6))
plt.plot(energy, dxy_up, label="dxy (up)", linestyle='-', linewidth=1.5)
plt.plot(energy, dxy_down, label="dxy (down)", linestyle='--', linewidth=1.5)

# Add labels, legend, and title
plt.xlabel("Energy (eV)")
plt.ylabel("Density of States (DOS)")
plt.title("Projected DOS for dxy Orbital (Gaussian Smearing)")
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.legend()
plt.grid(alpha=0.3)

# Display the plot
plt.tight_layout()
plt.show()
