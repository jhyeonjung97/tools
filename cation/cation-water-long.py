import os
import numpy as np
from ase.io import read, write
from ase import Atoms

path = '/Users/hailey/Desktop/1_cation/3_PES/1_cation_water/0_'
dirs = ['1_Li', '2_Na', '3_K', '4_Rb', '5_Cs']
for i, dir in enumerate(dirs):
    json_cation = os.path.join(path, dir, 'restart.json')
    atoms = read(json_cation)
    l1 = atoms.cell.lengths()[0]
    l2 = atoms.cell.lengths()[1]
    l3 = atoms.cell.lengths()[2]
    a1 = atoms.cell.angles()[0]
    a2 = atoms.cell.angles()[1]
    a3 = atoms.cell.angles()[2]
    if i == 0 or i == 1:
        atoms.positions[:, 2] -= l3/6
        atoms.wrap()
    atoms.positions[:, 2] += l3/3
    atoms.set_cell([l1, l2, l3/3*5, a1, a2, a3], scale_atoms=False)

    sorted_indices = np.argsort(atoms.positions[:, 2])
    top36_indices = sorted_indices[-36:]
    new_positions = atoms.positions[top36_indices].copy()
    new_positions[:, 2] += l3/3
    new_symbols = [atoms.symbols[i] for i in top36_indices]
    new_atoms = Atoms(symbols=new_symbols, positions=new_positions, cell=atoms.cell, pbc=atoms.pbc)
    atoms.extend(new_atoms)
    bottom36_indices = sorted_indices[:36]
    new_positions = atoms.positions[bottom36_indices].copy()
    new_positions[:, 2] -= l3/3
    new_symbols = [atoms.symbols[i] for i in bottom36_indices]
    new_atoms = Atoms(symbols=new_symbols, positions=new_positions, cell=atoms.cell, pbc=atoms.pbc)
    atoms.extend(new_atoms)

    write(f'restart_{i+1}.json', atoms)    

