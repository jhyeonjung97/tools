from typing import Any
from ase.io import read, write
from ase import Atoms
import os

for file in os.listdir('./'):
    if file.endswith('.json'):
        filename = file.split('.')[0]
        atoms = read(f'{filename}.json')
        min_z = atoms.positions[:,2].min()
        # max_z = atoms.positions[:,2].max()
        # height = max_z - min_z + 15.0
        atoms.positions[:, 2] -= min_z - 0.1

        l1 = atoms.cell.lengths()[0]
        l2 = atoms.cell.lengths()[1]
        a1 = atoms.cell.angles()[0]
        a2 = atoms.cell.angles()[1]
        a3 = atoms.cell.angles()[2]
        atoms.cell = (l1, l2, 20, a1, a2, a3)

        # atoms.positions[:, 1] -= 0.009
        # atoms.wrap()

        # atoms.center(axis=2)

        # symbols = atoms.get_chemical_symbols()
        # ir_indices = [i for i, sym in enumerate[Any](symbols) if sym == 'Ir']
        # fe_indices = [i for i, sym in enumerate[Any](symbols) if sym == 'Fe']
        # o_indices = [i for i, sym in enumerate[Any](symbols) if sym == 'O']
        # h_indices = [i for i, sym in enumerate[Any](symbols) if sym == 'H']
        # new_indices = ir_indices + fe_indices + o_indices + h_indices
        # atoms = atoms[new_indices]

        atoms.write(f'{filename}.cif', format='cif')