import os
import time
import shutil
import subprocess
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.vasp import Vasp

root = os.getcwd()
basename = os.path.basename(root)
dirs = ['0_H', '1_Li', '2_Na', '3_K', '4_Rb', '5_Cs']

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
    ncore=4,
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
    nelm=600,
    setups='recommended',
    # lasph=True,
    # laechg=True,
    # isym=0,
    lorbit=11,
    # lwave=False,
    )

for dir in dirs:

    os.chdir(os.path.join(root, dir, '00_'))
    atoms = read('restart.json')

    atoms.calc = calc
    atoms.get_potential_energy()

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