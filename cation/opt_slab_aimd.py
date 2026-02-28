import os
import time
import subprocess
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'opt_slab_aimd'
start_time = time.time()

if os.path.exists('restart.json'):
    atoms = read('restart.json')
else:
    atoms = read('start.traj')

atoms.calc = vasp_calculator.Vasp(
    istart=1,
    icharg=1,
    encut=400,
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
    ibrion=0,
    isif=2,
    ispin=1,
    ediffg=-0.02,
    ediff=1e-4,
    nsw=3000,
    potim=0.5,
    tebeg=300,
    teend=300,
    smass=0,
    setups='recommended',
    lasph=True,
    # laechg=True,
    # isym=0,
    # lorbit=11,
    # lwave=False,
    # lsol=True,
    )

eng = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
subprocess.call(f'python ~/bin/get_restart3', shell=True)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("="*40 + "\n")

print(f"Execution time logged in 'time.log'.")