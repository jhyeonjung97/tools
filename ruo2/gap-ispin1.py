from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin
import numpy as np
import matplotlib.pyplot as plt

# Load vasprun.xml file
vasprun = Vasprun("vasprun.xml", parse_dos=True)
complete_dos = vasprun.complete_dos

# Get DOS for non-spin-polarized calculation
densities = complete_dos.densities

# Extract energies and densities (only one spin channel for ISPIN=1)
energies = complete_dos.energies - complete_dos.efermi
total_dos = densities[Spin.up]  # For ISPIN=1, use only spin-up channel

# Determine the DOS grid resolution
dos_grid_resolution = energies[1] - energies[0]
# print(dos_grid_resolution)

# Find the band gap
conduction_band_min = None
valence_band_max = None

for i, energy in enumerate(energies):
    if total_dos[i] > 0 and energy < 0:
        valence_band_max = energy
    elif total_dos[i] > 0 and energy > 0 and conduction_band_min is None:
        conduction_band_min = energy
# print(valence_band_max, conduction_band_min)

if valence_band_max is not None and conduction_band_min is not None:
    band_gap = conduction_band_min - valence_band_max
    if band_gap <= dos_grid_resolution + 0.001:
        message = "Band Gap: 0.000 eV (metallic)"
    else:
        message = f"Band Gap: {band_gap:.3f} eV"
else:
    message = "No band gap found."
print(f'{message}')

# Save the message to gap.txt
with open("BANDGAP", "w") as file:
    file.write(message)
