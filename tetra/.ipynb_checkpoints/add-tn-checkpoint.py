import os
import shutil
import subprocess
import numpy as np
from ase.io import read, write

# Define root directory and coordination type
root = "/pscratch/sd/j/jiuy97/8_V_slab"
coord = '3_SquarePlanar_TN'

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

            # Get bond vector from atom 19 to 9
            bond_vector = atoms.positions[9] - atoms.positions[19]

            # 1) Add oxygen to atom 9 using bond_vector
            pos_9 = atoms.positions[9] + bond_vector
            atoms_o1 = atoms.copy()
            atoms_o1.append('O')
            atoms_o1.positions[-1] = pos_9
            write(os.path.join(path, 'o1.json'), atoms_o1)

            # Add oxygen to atom 11 using flipped y in bond_vector
            flipped_vector = bond_vector * np.array([1.0, -1.0, 1.0])
            pos_11 = atoms.positions[11] + flipped_vector
            atoms_o2 = atoms.copy()
            atoms_o2.append('O')
            atoms_o2.positions[-1] = pos_11
            write(os.path.join(path, 'o2.json'), atoms_o2)

            # 2) Add H on each oxygen in o1 and o2
            atoms_oh1 = atoms_o1.copy()
            h1_pos = pos_9 + np.array([0.0, 1.0, 0.0])
            atoms_oh1.append('H')
            atoms_oh1.positions[-1] = h1_pos
            write(os.path.join(path, 'oh1.json'), atoms_oh1)

            atoms_oh2 = atoms_o2.copy()
            h2_pos = pos_11 + np.array([0.0, -1.0, 0.0])
            atoms_oh2.append('H')
            atoms_oh2.positions[-1] = h2_pos
            write(os.path.join(path, 'oh2.json'), atoms_oh2)

            # 3) Make folders o1, o2, oh1, oh2 and move each as restart.json
            for name in ['o1', 'o2', 'oh1', 'oh2']:
                subdir = os.path.join(path, name)
                os.makedirs(subdir, exist_ok=True)
                shutil.copy(os.path.join(path, f'{name}.json'), os.path.join(subdir, 'restart.json'))

            # 4–5) Copy submit.sh and modify job names
            submit_path = os.path.join(path, 'submit.sh')
            if os.path.exists(submit_path):
                for suffix in ['o1', 'o2', 'oh1', 'oh2']:
                    subdir = os.path.join(path, suffix)
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
            subprocess.run(['sbatch', 'submit.sh'], cwd=folder, check=True)
            print(f"Submitted job in {folder}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit job in {folder}: {e}")

if __name__ == "__main__":
    main()
