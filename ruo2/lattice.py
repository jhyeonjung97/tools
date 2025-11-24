import os
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

parent_dir = '/pscratch/sd/j/jiuy97/3_RuO2/4_high_valence/1_M-RuO2/4_Re'
restart_path = parent_dir + '/restart.json'
wavecar_path = parent_dir + '/WAVECAR'
chgcar_path = parent_dir + '/CHGCAR'

for a in [6.43, 6.44]:
    for b in [6.44, 6.45, 6.46, 6.47, 6.48]:
        for c in [6.27]:
            folder_name = f'rutile_{a:.2f}_{b:.2f}_{c:.2f}'
            folder_path = os.path.join(parent_dir, 'others', folder_name)
            done_path = os.path.join(parent_dir, folder_name, 'DONE')
            if os.path.exists(done_path):
                continue
            os.makedirs(folder_path, exist_ok=True)
            # shutil.copy(restart_path, folder_path)
            # shutil.copy(wavecar_path, folder_path)
            # shutil.copy(chgcar_path, folder_path)
            os.chdir(folder_path)
            # atoms = read('restart.json')
            # atoms.set_cell([a+i, b+j, c+k])
            # atoms.write('restart.json')
            subprocess.call(f'sh ~/bin/verve/jobname.sh {folder_name}', shell=True)
            os.chdir('../')