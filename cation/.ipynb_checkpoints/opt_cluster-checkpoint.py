import os
import time
import subprocess
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'opt_cluster'
start_time = time.time()

if os.path.exists('restart.json'):
    atoms = read('restart.json')
else:
    atoms = read('start.traj')

atoms.calc = vasp_calculator.Vasp(
    encut=400,
    xc='PBE',
    gga='PE',
    kpts=(1,1,1),
    kpar=1,
    npar=1,
    gamma=True,
    ismear=0,
    sigma=0.05,
    algo='fast',
    lreal='auto',
    ibrion=2,
    isif=2,
    ispin=1,
    ediffg=-0.03,
    ediff=1e-4,
    nsw=800,
    setups='recommended',
    laechg=True,
    isym=0,
    # lsol=True,
    )

eng = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_{name}.json', shell=True)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Calculation Directory: {current_dir}\n")
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("="*40 + "\n")

print(f"Execution time logged in 'time.log'.")