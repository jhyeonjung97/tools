import subprocess
import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('restart.json')

atoms.calc = vasp_calculator.Vasp(
    amix=0.1, 
    amix_mag=0.1, 
    bmix=1e-05, 
    bmix_mag=1e-05, 
    encut=600, 
    sigma=0.05, 
    ediff=1e-05, 
    ediffg=-0.02, 
    algo='normal', 
    gga='PE', 
    prec='Normal', 
    ibrion=2, 
    isif=3, 
    ismear=0, 
    ispin=2, 
    istart=1, 
    # kpar=10,
    ldauprint=2, 
    ldautype=2, 
    lmaxmix=6, 
    lorbit=11, 
    nelm=250, 
    npar=6, 
    nsw=99, 
    # nupdown=0, 
    inimix=0, 
    laechg=True, 
    lasph=True, 
    ldau=True, 
    lvtot=False, 
    lreal=False, 
    ldau_luj={
        'Mn': {'L': 2, 'U': 2.75, 'J': 0.0},
        'O': {'L': -1, 'U': 0.0, 'J': 0.0},
        'H': {'L': -1, 'U': 0.0, 'J': 0.0}
        },
    xc='PBE',
    pp='PBE', 
    setups={'base': 'minimal'}, 
    txt='-', 
    kpts=(18, 5, 5), 
    gamma=True, 
    reciprocal=False, 
    ignore_constraints=False
    )

atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)