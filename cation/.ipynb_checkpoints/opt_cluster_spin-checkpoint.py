import subprocess
from os import path
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory

name = 'opt_cluster'

if path.exists('restart.json'):
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
    ispin=2,
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