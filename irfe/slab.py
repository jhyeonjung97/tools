from ase.io import read, write
from ase import Atoms
import os

for file in os.listdir('./'):
    if file.endswith('.json'):
        filename = file.split('.')[0]
        atoms = read(file)

        min_z = atoms.positions[:,2].min()
        max_z = atoms.positions[:,2].max()
        height = max_z - min_z + 15.0

        l1 = atoms.cell.lengths()[0]
        l2 = atoms.cell.lengths()[1]
        a1 = atoms.cell.angles()[0]
        a2 = atoms.cell.angles()[1]
        a3 = atoms.cell.angles()[2]
        atoms.cell = (l1, l2, height, a1, a2, a3)

        symbols = atoms.get_chemical_symbols()
        ir_indices = [i for i, sym in enumerate(symbols) if sym == 'Ir']
        o_indices = [i for i, sym in enumerate(symbols) if sym == 'O']
        other_indices = [i for i, sym in enumerate(symbols) if sym != 'Ir' and sym != 'O']
        new_indices = ir_indices + o_indices + other_indices
        atoms = atoms[new_indices]

        atoms.write(f'{filename}.vasp')