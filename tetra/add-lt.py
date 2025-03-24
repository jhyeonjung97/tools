import os
import shutil
import subprocess
import numpy as np
from ase.io import read, write
from ase.build import add_adsorbate

# Define root directory and coordination type
root = "/pscratch/sd/j/jiuy97/8_V_slab"
coord = '7_Pyramidal_LT'

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
            path = os.path.join(root, coord, row, f'{numb}_{metal}')
            atoms_path = os.path.join(path, 'final_with_calculator.json')
            
            if not os.path.exists(atoms_path):
                continue
    
            atoms = read(atoms_path)
    
            # 1) Add one oxygen atom 2 Å above the topmost atom in z
            top_z = atoms.positions[:, 2].max()
            top_index = np.argmax(atoms.positions[:, 2])
            top_position = atoms.positions[top_index]
            o_position = top_position + np.array([0, 0, 2.0])
            atoms_o = atoms.copy()
            atoms_o.append('O')
            atoms_o.positions[-1] = o_position
            write(os.path.join(path, 'o.json'), atoms_o)
    
            # 2) Add one hydrogen atom on the oxygen atom
            atoms_oh = atoms_o.copy()
            oh_vector = np.array([0.0, 0.8, 0.6])
            h_position = o_position + oh_vector
            atoms_oh.append('H')
            atoms_oh.positions[-1] = h_position
            write(os.path.join(path, 'oh.json'), atoms_oh)
    
            # 3) Create directories o/ and oh/
            o_dir = os.path.join(path, 'o')
            oh_dir = os.path.join(path, 'oh')
            os.makedirs(o_dir, exist_ok=True)
            os.makedirs(oh_dir, exist_ok=True)
    
            # 4) Put o.json and oh.json into respective directories as restart.json
            shutil.copy(os.path.join(path, 'o.json'), os.path.join(o_dir, 'restart.json'))
            shutil.copy(os.path.join(path, 'oh.json'), os.path.join(oh_dir, 'restart.json'))
    
            # 5) Copy submit.sh to both directories
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
            subprocess.run(['sbatch', 'submit.sh'], cwd=folder, check=True)
            print(f"Submitted job in {folder}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit job in {folder}: {e}")
            
if __name__ == "__main__":
    main()