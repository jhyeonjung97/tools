import os
import shutil
import subprocess
import numpy as np
import socket
from ase.io import read, write

# Define root directory and coordination type
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/8_V_slab'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/8_V_slab'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/8_V_slab'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")

coord = '4_SquarePlanar_PD'

rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

def main():
    for row in rows:
        metals = rows[row]
        for m, metal in enumerate(metals):
            numb = str(m).zfill(2)
            path = os.path.join(root, coord, row, f'{numb}_{metal}', 'clean')
            atoms_path = os.path.join(path, 'restart.json')
            
            if not os.path.exists(atoms_path) or os.path.exists(os.path.join(path, 'unmatched')):
                continue

            atoms = read(atoms_path)
            
            # Get positions of all atoms
            positions = atoms.positions
            symbols = atoms.get_chemical_symbols()
            
            # Separate oxygen and metal atoms
            o_indices = [i for i, sym in enumerate(symbols) if sym == 'O']
            metal_indices = [i for i, sym in enumerate(symbols) if sym != 'O']
            
            # Sort by z-coordinate
            o_positions = positions[o_indices]
            metal_positions = positions[metal_indices]
            
            o_z_sorted = o_positions[o_positions[:, 2].argsort()]
            metal_z_sorted = metal_positions[metal_positions[:, 2].argsort()]
            
            # Get lowest two metal atoms
            lowest_metals = metal_z_sorted[:2]
            # Find the metal atom with larger y-coordinate
            metal_with_larger_y = lowest_metals[np.argmax(lowest_metals[:, 1])]
            
            # Get third lowest oxygen
            third_lowest_o = o_z_sorted[2]
            
            # Calculate o1_vector (oxygen - metal)
            o1_vector = third_lowest_o - metal_with_larger_y
            bond_length = np.linalg.norm(o1_vector)
            
            # Get highest two metal atoms
            highest_metals = metal_z_sorted[-2:]
            # Calculate midpoint using cell's x-length
            x_midpoint = atoms.cell[0, 0] / 2
            # Find the metal atom that is further from the midpoint
            x_distances = np.abs(highest_metals[:, 0] - x_midpoint)
            further_metal = highest_metals[np.argmax(x_distances)]
            closer_metal = highest_metals[np.argmin(x_distances)]
            
            # Add first oxygen using o1_vector to the further metal
            pos_o1 = further_metal + o1_vector
            atoms_o1 = atoms.copy()
            atoms_o1.append('O')
            atoms_o1.positions[-1] = pos_o1
            write(os.path.join(path, 'o1.json'), atoms_o1)
            
            # Calculate o2_vector by flipping o1_vector in y-direction
            a = atoms.cell.cellpar()[0]
            b = atoms.cell.cellpar()[1]
            o2_vector = np.array([a/4, b/4, o1_vector[2]-0.5])
            
            # Add second oxygen using o2_vector to the closer metal
            pos_o2 = further_metal + o2_vector
            atoms_o2 = atoms.copy()
            atoms_o2.append('O')
            atoms_o2.positions[-1] = pos_o2
            write(os.path.join(path, 'o2.json'), atoms_o2)
            
            # Add hydrogen atoms in y-direction
            atoms_oh1 = atoms_o1.copy()
            h1_pos = pos_o1 + np.array([1.0, 0.0, 0.0])
            atoms_oh1.append('H')
            atoms_oh1.positions[-1] = h1_pos
            write(os.path.join(path, 'oh1.json'), atoms_oh1)
            
            atoms_oh2 = atoms_o2.copy()
            h2_pos = pos_o2 + np.array([0.0, 0.0, 1.0])
            atoms_oh2.append('H')
            atoms_oh2.positions[-1] = h2_pos
            write(os.path.join(path, 'oh2.json'), atoms_oh2)
            
            # Create directories and move files
            for name in ['o1', 'o2', 'oh1', 'oh2']:
                subdir = os.path.join(root, coord, row, f'{numb}_{metal}', name)
                os.makedirs(subdir, exist_ok=True)
                shutil.copy(os.path.join(path, f'{name}.json'), os.path.join(subdir, 'restart.json'))
            
            # Copy submit.sh and modify job names
            submit_path = os.path.join(path, 'submit.sh')
            if os.path.exists(submit_path):
                for suffix in ['o1', 'o2', 'oh1', 'oh2']:
                    subdir = os.path.join(root, coord, row, f'{numb}_{metal}', suffix)
                    shutil.copy(submit_path, os.path.join(subdir, 'submit.sh'))
                    modify_job_name(os.path.join(subdir, 'submit.sh'), suffix)
                    submit_job(subdir)

# Update the job name in submit.sh files
def modify_job_name(file_path, suffix):
    if not os.path.exists(file_path):
        return
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        for line in lines:
            if line.strip().startswith('#SBATCH -J'):
                line = line.rstrip() + f"{suffix}\n"
            f.write(line)

# Submit the job in each folder
def submit_job(folder):
    submit_script = os.path.join(folder, 'submit.sh')
    if os.path.exists(submit_script):
        try:
            # Read and modify submit.sh
            with open(submit_script, 'r') as f:
                lines = f.readlines()
            
            with open(submit_script, 'w') as f:
                for line in lines:
                    if line.startswith('#SBATCH -N'):
                        line = '#SBATCH -N 2\n'
                    elif line.startswith('#SBATCH -q'):
                        line = '#SBATCH -q regular\n'
                    elif line.startswith('#SBATCH -t'):
                        line = '#SBATCH -t 04:00:00\n'
                    elif line.startswith('#SBATCH -G'):
                        line = '#SBATCH -G 8\n'
                    elif 'run_vasp_gpu' in line:
                        line = line.replace('run_vasp_gpu.py', 'run_vasp_gpu2.py')
                        line = line.replace('run_vasp_gpu1.py', 'run_vasp_gpu2.py')
                        line = line.replace('run_vasp_gpu3.py', 'run_vasp_gpu2.py')
                        line = line.replace('run_vasp_gpu4.py', 'run_vasp_gpu2.py')
                        line = line.replace('run_vasp_gpu8.py', 'run_vasp_gpu2.py')
                        line = line.replace('run_vasp_gpu16.py', 'run_vasp_gpu2.py')
                    elif 'cpu' in line:
                        line = line.replace('cpu', 'gpu')
                    elif 'opt_slab2_afm.py' in line:
                        line = line.replace('opt_slab2_afm.py', 'opt_ads_gpu.py')
                    elif 'opt_slab2_fm.py' in line:
                        line = line.replace('opt_slab2_fm.py', 'opt_ads_gpu.py')
                    f.write(line)
            
            # Check if run_vasp_gpu2.py is present in the modified file
            with open(submit_script, 'r') as f:
                content = f.read()
                if 'run_vasp_gpu2.py' not in content:
                    print(f"Warning: run_vasp_gpu2.py not found in {submit_script}")
                    return
            
            # subprocess.run(['sbatch', 'submit.sh'], cwd=folder, check=True)
            # print(f"Submitted job in {folder}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit job in {folder}: {e}")

if __name__ == "__main__":
    main()
