from ase.io import read, write
from ase.geometry import cellpar_to_cell

atoms = read('POSCAR')
a = atoms.cell.lengths()[0]  # Use the first cell length
new_cell_parameters = [a, a, a, 100/3, 100/3, 100/3]
new_cell = cellpar_to_cell(new_cell_parameters)
atoms.set_cell(new_cell, scale_atoms=True)
write('start.traj', atoms)
