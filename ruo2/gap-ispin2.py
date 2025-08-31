from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin
import numpy as np
import matplotlib.pyplot as plt

# Load vasprun.xml file
vasprun = Vasprun("vasprun.xml", parse_dos=True)
complete_dos = vasprun.complete_dos

# Get spin-polarized DOS
densities = complete_dos.densities

# Extract energies and densities for spin-up and spin-down
energies = complete_dos.energies - complete_dos.efermi
spin_up_dos = densities[Spin.up]
spin_down_dos = densities[Spin.down]

# Determine the DOS grid resolution
dos_grid_resolution = energies[1] - energies[0]
# print(dos_grid_resolution)

# Find the band gap
conduction_band_min = None
valence_band_max = None

for i, energy in enumerate(energies):
    if (spin_up_dos[i] > 0 or spin_down_dos[i] > 0) and energy < 0:
        valence_band_max = energy
    elif (spin_up_dos[i] > 0 or spin_down_dos[i] > 0) and energy > 0 and conduction_band_min is None:
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
