import os
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

ldau_luj = {
    'Ru': {'L': 2, 'U': 2.5, 'J': 0.0},
    'Re': {'L': -1, 'U': 0.0, 'J': 0.0},
    'O': {'L': -1, 'U': 0.0, 'J': 0.0},
    'H': {'L': -1, 'U': 0.0, 'J': 0.0}
    }

calc = vasp_calculator.Vasp(
    #inimix=0,
    #amix=0.1,
    #amix_mag=0.1,
    #bmix=0.0001,
    #bmix_mag=0.0001,
    encut=500,
    sigma=0.05,
    ediff=1e-05,
    ediffg=-0.02,
    algo='Fast',
    gga='PE',
    prec='Normal',
    ibrion=2,
    isif=2,
    ismear=0,
    ispin=2,
    istart=1,
    ldau=True,
    lmaxmix=4,
    ldau_luj=ldau_luj,
    ldautype=2,
    ldauprint=2,
    lorbit=11,
    nelm=250,
    npar=6,
    nsw=199,
    laechg=True,
    lasph=True,
    lvtot=False,
    lreal='Auto',
    kpts=[5, 5, 5],
    gamma=True,
    )

parent_dir = '/pscratch/sd/j/jiuy97/3_RuO2/4_high_valence/1_M-RuO2/0_Ru'
restart_path = parent_dir + '/restart.json'
wavecar_path = parent_dir + '/WAVECAR'

a = 6.43
b = 6.46
c = 6.24

lattice_change = [-0.1, 0.0, +0.1]
for i in lattice_change:
    for j in lattice_change:
        for k in lattice_change:
            folder_name = f'rutile_{a+i:.2f}_{b+j:.2f}_{c+k:.2f}'
            folder_path = os.path.join(parent_dir, folder_name)
            os.makedirs(folder_path)
            shutil.copy(restart_path, folder_path)
            shutil.copy(wavecar_path, folder_path)
            os.chdir(folder_path)
            atoms = read('restart.json')
            atoms.set_cell([a+i, b+j, c+k])
            atoms.calc = calc
            atoms.write('start.traj')
            atoms.get_potential_energy()
            print('Calculation Complete, storing the run + calculator to traj file')
            subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)
            Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
            subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
            os.chdir('../')