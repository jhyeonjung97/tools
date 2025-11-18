import os
import mendeleev as md
from ase.io import read

for i in list(range(21, 31)) + list(range(39, 49)) + [57] + list(range(72, 81)):
    element = md.element(i).symbol
    dir = f'{i}_{element}'
    os.makedirs(dir, exist_ok=True)
    atoms = read(f'restart.json')
    for i in range(0, 8):
        atoms[i].symbol = element
    atoms.write(f'{dir}/restart.json')