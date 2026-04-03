import os
import time
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.vasp import Vasp

root = os.getcwd()
basename = os.path.basename(root)

name = 'pes_slab'
start_time = time.time()

calc = Vasp(
    istart=1,
    icharg=1,
    inimix=0,
    amix=0.1,
    bmix=0.0001,
    amix_mag=0.1,
    bmix_mag=0.0001,
    encut=500,
    xc='PBE',
    gga='RP',
    ivdw=12,
    kpts=(3, 3, 1),
    kpar=5,
    ncore=6,
    gamma=True,
    ismear=0,
    sigma=0.05,
    algo='normal',
    lreal='auto',
    ibrion=-1,
    isif=2,
    ispin=2,
    ediffg=-0.02,
    ediff=1e-5,
    nsw=0,
    nelm=250,
    setups='recommended',
    # lasph=True,
    # laechg=True,
    # isym=0,
    lorbit=11,
    # lwave=False,
    )

for i in range(30):
    j = i + 1

    old_folder_name = f'{i:02d}_'
    old_folder_path = os.path.join(root, old_folder_name)

    new_folder_name = f'{j:02d}_'
    new_folder_path = os.path.join(root, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)

    shutil.copy2(os.path.join(old_folder_path, 'restart.json'), new_folder_path)  
    shutil.copy2(os.path.join(old_folder_path, 'WAVECAR'), new_folder_path)
    shutil.copy2(os.path.join(old_folder_path, 'CHGCAR'), new_folder_path)

    os.chdir(new_folder_path)
    atoms = read('restart.json')
    z_coords = atoms.get_positions()[:, 2]
    min_z_index = z_coords.argmin()
    atoms[min_z_index].z -= 0.1
    write('restart.json', atoms)

    atoms.calc = calc
    atoms.get_potential_energy()
    atoms.get_forces()

    Trajectory(f'final_{name}.traj', 'w').write(atoms)
    subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
    subprocess.call(f'/sdf/home/j/jiuy97/miniforge3/bin/python ~/bin/get_restart3', shell=True)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("="*40 + "\n")

print(f"Execution time logged in 'time.log'.")