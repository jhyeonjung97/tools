import os
import shutil
import subprocess
import numpy as np
import socket
from ase.io import read, write
from ase.build import add_adsorbate

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

coord = '9_Tetrahedral_ZB'

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
            
            # Get second highest oxygen and third highest metal
            second_highest_o = o_z_sorted[-2]
            third_highest_metal = metal_z_sorted[-3]
            
            # Calculate o_vector
            o_vector = second_highest_o - third_highest_metal
            
            # Get highest metal position and add new oxygen
            highest_metal = metal_z_sorted[-1]
            new_o_position = highest_metal + o_vector
            
            # Create new structure with additional oxygen
            atoms_o = atoms.copy()
            atoms_o.append('O')
            atoms_o.positions[-1] = new_o_position
            write(os.path.join(path, 'o.json'), atoms_o)
            
            # Add hydrogen atom 0.9 Angstroms in x-direction from the new oxygen
            atoms_oh = atoms_o.copy()
            h_position = new_o_position + np.array([0.9, 0.0, 0.0])
            atoms_oh.append('H')
            atoms_oh.positions[-1] = h_position
            write(os.path.join(path, 'oh.json'), atoms_oh)
    
            # Create directories o/ and oh/
            o_dir = os.path.join(root, coord, row, f'{numb}_{metal}', 'o')
            oh_dir = os.path.join(root, coord, row, f'{numb}_{metal}', 'oh')
            os.makedirs(o_dir, exist_ok=True)
            os.makedirs(oh_dir, exist_ok=True)
    
            # Put o.json and oh.json into respective directories as restart.json
            shutil.copy(os.path.join(path, 'o.json'), os.path.join(o_dir, 'restart.json'))
            shutil.copy(os.path.join(path, 'oh.json'), os.path.join(oh_dir, 'restart.json'))
    
            # Copy submit.sh to both directories
            submit_path = os.path.join(path, 'submit.sh')
            if os.path.exists(submit_path):
                shutil.copy(submit_path, o_dir)
                shutil.copy(submit_path, oh_dir)
                modify_job_name(os.path.join(o_dir, 'submit.sh'), 'o')
                modify_job_name(os.path.join(oh_dir, 'submit.sh'), 'oh')
                submit_job(o_dir)
                submit_job(oh_dir)

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
                    elif 'WZ' in line:
                        line = line.replace('ZB', 'ZZ')
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