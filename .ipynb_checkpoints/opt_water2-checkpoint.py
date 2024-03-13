import sys
import json
import subprocess
import numpy as np
from os import path
from mendeleev import element
from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'opt_water2'

if path.exists('restart.json'):
    atoms = read('restart.json')
else:
    atoms = read('start.traj')
    i = 1

atoms.calc = vasp_calculator.Vasp(
                    encut=400,
                    xc='PBE',
                    gga='PE',
                    kpts=(1,1,1),
                    kpar=8,
                    npar=1,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    algo='fast',
                    lreal='auto',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.03,
                    ediff=1e-4,
                    nsw=800,
                    setups='recommended',
                    laechg=True,
                    isym=0
                    )

eng = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'restart.traj','w').write(atoms)
subprocess.call(f'cp restart.traj final_{name}.traj', shell=True)
subprocess.call(f'ase convert -f final_{name}.traj final_{name}.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)