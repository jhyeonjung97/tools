import os
from ase.io import read, write

for file in os.listdir('.'):
    if file.endswith('.json'):
        atoms = read(file)
        write(file.replace('.json', '.vasp'), atoms)