import os
import glob
import re
from ase.io import read, write

for i in range(8):
    contcar_path=f'/pscratch/sd/j/jiuy97/6_MNC/empty/{i}_/CONTCAR'
    atoms = read(contcar_path)
    vasp_output_file = f'/pscratch/sd/j/jiuy97/6_MNC/figures/vasp/vac{i}.vasp'
    write(vasp_output_file, atoms)
    print(f"Saved as vac{i}.vasp")
    png_file = f'/pscratch/sd/j/jiuy97/6_MNC/figures/vasp/vac{i}.png'
    write(png_file, atoms, rotation='-90x, -90y, 0z', show_unit_cell=0)
    print(f"Saved as vac{i}.png")
