from ase.io import read, write
import glob
import os

for i in range(7): # 0~6
    path_pattern = f'/scratch/x2755a09/3_MNC/*d/*_*/*_*S/{i}_'
    matching_paths = glob.glob(path_pattern)
    for path in matching_paths:
        poscar_path = os.path.join(path, 'POSCAR')
        if os.path.exists(poscar_path):
            atoms = read(poscar_path)
            for atom in atoms:
                if atom.symbol == 'N':
                    atom.z = 10.0
                elif atom.symbol not in ['C', 'N', 'O', 'H']:
                    atom.z = 10.0 + i * 0.2
            output_path = os.path.join(path, 'start.traj')
            write(output_path, atoms)
        else:
            print(f'POSCAR file not found in {path}')