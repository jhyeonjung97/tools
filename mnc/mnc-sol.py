import sys
import json
import subprocess
import numpy as np
from os import path
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'mnc'
effective_length = 25

lmaxmix = 2
ldau_luj = {
    'Ti': {'L': 2, 'U': 3.00, 'J': 0.0},
    'V': {'L': 2, 'U': 3.25, 'J': 0.0},
    'Cr': {'L': 2, 'U': 3.5, 'J': 0.0},
    'Mn': {'L': 2, 'U': 3.75, 'J': 0.0},
    'Fe': {'L': 2, 'U': 4.3, 'J': 0.0},
    'Co': {'L': 2, 'U': 3.32, 'J': 0.0},
    'Ni': {'L': 2, 'U': 6.45, 'J': 0.0},
    'Cu': {'L': 2, 'U': 3.0, 'J': 0.0},
}

if path.exists('restart.json'):
    atoms = read('restart.json')
    with open('restart.json', 'r') as file:
        data = json.load(file)
    calculator_parameters = data['1']['calculator_parameters']
    if calculator_parameters:
        atoms.calc = vasp_calculator.Vasp(**calculator_parameters)
    else:
        for atom in atoms:
            if atom.symbol in ldau_luj:
                lmaxmix = 4
            else:
                ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
        atoms.calc = vasp_calculator.Vasp(
            istart=0,
            encut=500,
            gga='PE',
            ivdw=12,
            kpts=(5, 5, 1),
            kpar=8,
            npar=1,
            gamma=True,
            ismear=0,
            sigma=0.05,
            algo='Normal',
            ibrion=2,
            isif=2,
            ediffg=-0.02,
            ediff=1e-5,
            prec='Normal',
            nsw=200,
            lvhar=True,
            lvtot=False,
            ispin=2,
            setups={'base': 'recommended', 'W': '_sv'},
            ldau=True,
            ldautype=2,
            ldau_luj=ldau_luj,
            ldauprint=2,
            lmaxmix=lmaxmix,
            lasph=True,
            laechg=True,
            lreal='Auto',
            nedos=3000,
            lorbit=11,
            lsol=True
        )
else:
    raise ValueError('restart.json not found')

e = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

traj_filename = f'final_{name}.traj'
Trajectory(traj_filename, 'w').write(atoms)
subprocess.call(f'ase convert -f {traj_filename} final_with_calculator.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)