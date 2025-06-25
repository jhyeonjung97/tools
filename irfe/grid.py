from ase.io import read
import numpy as np

atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()
for i in range(5):
    scale_x = 1.0 + (i) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x    
    for j in range(5):
        scale_y = 1.0 + (j) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y        
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y        
        atoms.write(f'a{i}{j}.vasp')

atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()
for i in range(5):
    scale_x = 1.0 - (i) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x    
    for j in range(5):
        scale_y = 1.0 + (j) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y        
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y        
        atoms.write(f'b{i}{j}.vasp')

atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()
for i in range(5):
    scale_x = 1.0 - (i) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x    
    for j in range(5):
        scale_y = 1.0 - (j) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y        
        atoms.write(f'c{i}{j}.vasp')

atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()
for i in range(5):
    scale_x = 1.0 + (i) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x    
    for j in range(5):
        scale_y = 1.0 - (j) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y        
        atoms.write(f'd{i}{j}.vasp')

# mkdir 00_ 01_ 02_ 03_ 04_
# mkdir 10_ 11_ 12_ 13_ 14_
# mkdir 20_ 21_ 22_ 23_ 24_
# mkdir 30_ 31_ 32_ 33_ 34_
# mkdir 40_ 41_ 42_ 43_ 44_