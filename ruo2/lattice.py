import os
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

parent_dir = '/pscratch/sd/j/jiuy97/3_RuO2/4_high_valence/1_M-RuO2/4_Re'
# parent_dir = '/pscratch/sd/j/jiuy97/3_RuO2/4_high_valence/1_M-RuO2/0_Ru'
restart_path = parent_dir + '/restart.json'
wavecar_path = parent_dir + '/WAVECAR'
chgcar_path = parent_dir + '/CHGCAR'
submit_path = parent_dir + '/submit.sh'
bulk_path = parent_dir + '/bulk.py'

a = 6.43
b = 6.48
c = 6.24

lattice_change = [-0.01, 0.00, +0.01]
for i in lattice_change:
    for j in lattice_change:
        for k in lattice_change:
            folder_name = f'rutile_{a+i:.2f}_{b+j:.2f}_{c+k:.2f}'
            folder_path = os.path.join(parent_dir, 'others', folder_name)
            done_path = os.path.join(parent_dir, folder_name, 'DONE')
            if os.path.exists(done_path):
                continue
            os.makedirs(folder_path, exist_ok=True)
            shutil.copy(restart_path, folder_path)
            shutil.copy(wavecar_path, folder_path)
            shutil.copy(chgcar_path, folder_path)
            shutil.copy(submit_path, folder_path)
            shutil.copy(bulk_path, folder_path)
            os.chdir(folder_path)
            atoms = read('restart.json')
            atoms.set_cell([a+i, b+j, c+k], scale_atoms=True)
            atoms.write('restart.json')
            subprocess.call(f'sh ~/bin/verve/jobname.sh {folder_name}', shell=True)
            os.chdir('../')

# for a in [6.46]:
#     for c in [6.21, 6.22, 6.24, 6.25]:
#         folder_name = f'rutile_{a:.2f}_{c:.2f}'
#         folder_path = os.path.join(parent_dir, 'others', folder_name)
#         done_path = os.path.join(parent_dir, folder_name, 'DONE')
#         if os.path.exists(done_path):
#             continue
#         os.makedirs(folder_path, exist_ok=True)
#         shutil.copy(restart_path, folder_path)
#         shutil.copy(wavecar_path, folder_path)
#         shutil.copy(chgcar_path, folder_path)
#         shutil.copy(submit_path, folder_path)
#         shutil.copy(bulk_path, folder_path)
#         os.chdir(folder_path)
#         atoms = read('restart.json')
#         atoms.set_cell([a, a, c], scale_atoms=True)
#         atoms.write('restart.json')
#         subprocess.call(f'sh ~/bin/verve/jobname.sh {folder_name}', shell=True)
#         os.chdir('../')