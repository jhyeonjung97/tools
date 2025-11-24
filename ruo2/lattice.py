import os
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

a = 6.43
b = 6.46
c = 6.24

parent_dir = '/pscratch/sd/j/jiuy97/3_RuO2/4_high_valence/1_M-RuO2/4_Re'
restart_path = parent_dir + '/restart.json'
wavecar_path = parent_dir + '/WAVECAR'
chgcar_path = parent_dir + '/CHGCAR'

lattice_change = [-0.02, -0.01, 0.00, +0.01, +0.02]
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
            os.chdir(folder_path)
            atoms = read('restart.json')
            atoms.set_cell([a+i, b+j, c+k])
            atoms.write('restart.json')
            os.chdir('../')