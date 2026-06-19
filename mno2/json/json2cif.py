from ase.io import read, write
import os

for file in os.listdir('./'):
    if file.endswith('.json'):
        atoms = read(file)
        write(file.replace('.json', '.cif'), atoms, format='cif')
