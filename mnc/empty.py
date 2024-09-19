from ase.io import read, write
from ase import Atoms

metal_index = None
atoms = read('../1_LS/relaxed/restart.json')

if atoms[-1].symbol not in ['C', 'N', 'O', 'H']:
    metal_index = atoms[-1].index

if metal_index is not None:
    del atoms[metal_index]
    write('restart.json', atoms)
