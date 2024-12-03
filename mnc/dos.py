import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Plot projected DOS without Gaussian smearing or Fermi level adjustment.")
parser.add_argument("--file", type=str, required=True, help="Path to the input DOS data file.")
parser.add_argument("--xrange", type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=(-8, 6),
                    help="Energy range for the x-axis (default: -8 to 6).")
args = parser.parse_args()

# Define column names
columns = [
    "energy", "s(up)", "s(down)", "p(up)", "p(down)",
    "dxy(up)", "dxy(down)", "dyz(up)", "dyz(down)",
    "dz2(up)", "dz2(down)", "dxz(up)", "dxz(down)",
    "dx2(up)", "dx2(down)"
]

# Read the DOS data
data = pd.read_csv(args.file, delim_whitespace=True, comment='#', names=columns)

# Extract energy data
energy = data["energy"]

# Define orbitals to plot
orbitals = ["dxy", "dyz", "dz2", "dxz", "dx2"]

# Create subplots
fig, axes = plt.subplots(len(orbitals), 1, figsize=(8, len(orbitals) * 2.5), sharex=True)

for i, orbital in enumerate(orbitals):
    orbital_up = data[f"{orbital}(up)"]
    orbital_down = data[f"{orbital}(down)"]

    # Plot on the respective subplot
    axes[i].plot(energy, orbital_up, label=f"{orbital} (up)", linestyle='-', linewidth=1.5)
    axes[i].plot(energy, orbital_down, label=f"{orbital} (down)", linestyle='--', linewidth=1.5)
    axes[i].fill_between(energy, 0, orbital_up, alpha=0.3)
    axes[i].fill_between(energy, 0, orbital_down, alpha=0.3)
    axes[i].axhline(0, color='black', linewidth=0.8, linestyle='--')
    axes[i].legend()
    axes[i].set_ylabel("DOS (arb. units)")
    axes[i].grid(alpha=0.3)

# Set shared x-axis labels
axes[-1].set_xlabel("Energy (eV)")
plt.xlim(args.xrange)

# Adjust layout
plt.tight_layout()
plt.show()
