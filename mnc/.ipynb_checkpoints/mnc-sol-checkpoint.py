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

# Define the name and effective length
name = 'mnc'
effective_length = 25

# Define LDA+U parameters for specific elements
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

# Check if the restart.json file exists
if path.exists('restart.json'):
    # Read the Atoms object from the restart file
    atoms = read('restart.json')

    try:
        # Open and read the JSON file to extract calculator parameters
        with open('restart.json', 'r') as file:
            data = json.load(file)
        calculator_parameters = data['1']['calculator_parameters']
        
        # Assign the calculator to the Atoms object
        atoms.calc = vasp_calculator.Vasp(**calculator_parameters)
    except Exception as e:
        print(f"Error reading calculator parameters: {e}")
        
        # Set default magnetic moments and calculator parameters if extraction fails
        magmoms = atoms.get_magnetic_moments()
        for atom in atoms:
            atom.magmom = magmoms[atom.index]
            if atom.symbol in ldau_luj:
                lmaxmix = 4
            else:
                ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
        
        # Set default VASP calculator parameters
        atoms.calc = vasp_calculator.Vasp(
            istart=0,
            encut=500,
            gga='PE',
            ivdw=12,
            kpts=(5, 5, 1),
            kpar=4,
            npar=16,
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

# Continue the calculation
e = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')

# Run post-processing scripts
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

# Save the final Atoms object with the calculator to a trajectory file
traj_filename = f'final_{name}.traj'
Trajectory(traj_filename, 'w').write(atoms)

# Convert the trajectory file to JSON format and copy necessary files
subprocess.call(f'ase convert -f {traj_filename} restart.json', shell=True)
subprocess.call(f'cp restart.json final_with_calculator.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)