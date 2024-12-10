import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Plot projected DOS for orbitals with unified y-axis and customizable y-title.")
parser.add_argument("--file", type=str, required=True, help="Path to the input DOS data file.")
parser.add_argument("--xrange", type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=(-8, 6),
                    help="Energy range for the x-axis (default: -8 to 6).")
parser.add_argument("--output", type=str, default=None, help="Output file name to save the plot.")
args = parser.parse_args()

# Define column names
# columns = [
#     "energy", "dxy(up)", "dxy(down)", "dyz(up)", "dyz(down)",
#     "dz2(up)", "dz2(down)", "dxz(up)", "dxz(down)", "dx2(up)", "dx2(down)"
# ]
columns = [
    "energy", "dx2(up)", "dx2(down)", "dyz(up)", "dyz(down)",
    "dz2(up)", "dz2(down)", "dxz(up)", "dxz(down)", "dxy(up)", "dxy(down)"
]
colors = ['blue', 'orange', 'green', 'red', 'purple']

# Read the DOS data
data = pd.read_csv(args.file, delim_whitespace=True, comment='#', names=columns)

# Extract energy data
energy = data["energy"]

# Define orbitals to plot
# orbitals = ["dxy", "dyz", "dz2", "dxz", "dx2"]
orbitals = ["dx2", "dxy", "dz2", "dxz", "dyz"]

# Find global y-axis range
y_min = float("inf")
y_max = float("-inf")

# Function to calculate d-band center
def calculate_d_band_center(energy, dos_up, dos_down):
    dos_total = dos_up - dos_down
    numerator = np.trapz(energy * dos_total, energy)
    denominator = np.trapz(dos_total, energy)
    return numerator / denominator

# Calculate d-band centers
d_band_centers = {}

for orbital in orbitals:
    dos_up = data[f"{orbital}(up)"]
    dos_down = data[f"{orbital}(down)"]
    d_band_centers[orbital] = calculate_d_band_center(energy, dos_up, dos_down)
    y_min = min(y_min, data[f"{orbital}(up)"].min(), data[f"{orbital}(down)"].min())
    y_max = max(y_max, data[f"{orbital}(up)"].max(), data[f"{orbital}(down)"].max())

# Save d-band centers to TSV
d_band_centers_df = pd.DataFrame(list(d_band_centers.items()), columns=["Orbital", "d_band_center"])
d_band_centers_df.to_csv(f"/pscratch/sd/j/jiuy97/6_MNC/figures/dos/{args.output}.tsv", sep="\t", index=False)
print(f"d-Band centers saved to {args.output}.tsv")
    
# Scale the limit symmetrically
ylimit = max(abs(y_min), abs(y_max)) * 1.2

# Create a single figure with stacked plots, no space between subplots
fig, axes = plt.subplots(len(orbitals), 1, figsize=(4, len(orbitals) * 1.5), sharex=True,
                         gridspec_kw={"hspace": 0})  # No vertical space between subplots

                
for i, orbital in enumerate(orbitals):
    orbital_up = data[f"{orbital}(up)"]
    orbital_down = data[f"{orbital}(down)"]

    # Plot on the respective subplot
    axes[i].plot(energy, orbital_up, label=f"{orbital}", color=colors[i], linestyle='-', linewidth=1.5)
    axes[i].plot(energy, orbital_down, label=None, color=colors[i], linestyle='-', linewidth=1.5)
    axes[i].fill_between(energy, 0, orbital_up, color=colors[i], alpha=0.3)
    axes[i].fill_between(energy, 0, orbital_down, color=colors[i], alpha=0.3)
    
    # Add horizontal and vertical dashed lines
    axes[i].axhline(0, color='black', linewidth=0.8, linestyle='--')  # Horizontal line at y=0
    axes[i].axvline(0, color='black', linewidth=0.8, linestyle='--')  # Vertical line at x=0
    
    # Set y-axis limits symmetrically
    axes[i].set_ylim(-ylimit, +ylimit)
    axes[i].yaxis.set_visible(False)  # Hide y-axis for other subplots
    if i != 4:
        axes[i].xaxis.set_visible(False)  # Hide y-axis for other subplots

    # Add legend
    axes[i].legend(loc="upper left", fontsize=10)

# Set shared x-axis labels
axes[-1].set_xlabel("Energy (eV)", fontsize=10)

# Set x-axis range
plt.xlim(args.xrange)

# Adjust layout
if args.output:
    plt.savefig(f"/pscratch/sd/j/jiuy97/6_MNC/figures/dos/{args.output}.png", dpi=300, bbox_inches='tight')  # Save with high resolution
    print(f"Plot saved as {args.output}.png")
else:
    plt.show()
