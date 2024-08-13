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

name = 'mnc-hs'

spin_states_plus_1 = {'Ti': 3, 'V': 4, 'Cr': 5, 'Mn': 4, 'Fe': 3, 'Co': 2,
                      'Zr': 3, 'Nb': 4, 'Mo': 5, 'Tc': 4, 'Ru': 3, 'Rh': 2,
                      'Hf': 3, 'Ta': 4, 'W': 5, 'Re': 4, 'Os': 3, 'Ir': 2
                      }
spin_states_plus_2 = {'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 4, 'Co': 3, 'Ni': 2,
                      'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 4, 'Rh': 3, 'Pd': 2,
                      'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 4, 'Ir': 3, 'Pt': 2
                      }
spin_states_plus_3 = {'V': 2, 'Cr': 3, 'Mn': 4, 'Fe': 5, 'Co': 4, 'Ni': 3, 'Cu': 2,
                      'Nb': 2, 'Mo': 3, 'Tc': 4, 'Ru': 5, 'Rh': 4, 'Pd': 3,
                      'Ta': 2, 'W': 3, 'Re': 4, 'Os': 5, 'Ir': 4, 'Pt': 3
                      }
spin_states_plus_4 = {'Cr': 2, 'Mn': 3, 'Fe': 4, 'Co': 5, 'Ni': 4, 'Cu': 3,
                      'Mo': 2, 'Tc': 3, 'Ru': 4, 'Rh': 5, 'Pd': 4,
                      'W': 2, 'Re': 3, 'Os': 4, 'Ir': 5, 'Pt': 4
                      }

ldau_luj = {'Ti': {'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr': {'L':2, 'U':3.50, 'J':0.0},
            'Mn': {'L':2, 'U':3.75, 'J':0.0},
            'Fe': {'L':2, 'U':4.30, 'J':0.0},
            'Co': {'L':2, 'U':3.32, 'J':0.0},
            'Ni': {'L':2, 'U':6.45, 'J':0.0},
            'Cu': {'L':2, 'U':3.00, 'J':0.0}
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
    amix_mag = 0.01
    bmix_mag = 0.00001
elif path.exists('start.traj'):
    atoms = read('start.traj')
    amix_mag = 0.05
    bmix_mag = 0.0001
    if atoms[-2].symbol == 'C' and atoms[-1].symbol == 'O':
        spin_states = spin_states_plus_2
    elif atoms[-2].symbol == 'O' and atoms[-1].symbol == 'H':
        spin_states = spin_states_plus_3
    elif atoms[-1].symbol == 'O':
        spin_states = spin_states_plus_4
    elif atoms[-1].symbol == 'H':
        spin_states = spin_states_plus_1
    else:
        spin_states = spin_states_plus_2
    for atom in atoms:
        if atom.symbol not in ['C', 'N', 'O', 'H']:
            if atom.symbol in spin_states:
                spin = spin_states.get(atom.symbol)
                atom.magmom = spin
            else:
                raise ValueError(f"Unexpected atom symbol '{atom.symbol}' found in start.traj")
else:
    raise ValueError('Neither restart.json nor start.traj file found')
        
lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    elif atom.symbol in spin_states:
        ldau_luj[atom.symbol] = {'L': 2, 'U': 0.0, 'J': 0.0}
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
        
atoms.calc = vasp_calculator.Vasp(
                    encut=500,
                    gga='PE',
                    ivdw=12,
                    kpts=(5,5,1),
                    kpar=8,
                    npar=1,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    inimix=0,
                    amix=0.05,
                    bmix=0.0001,
                    amix_mag=amix_mag,
                    bmix_mag=bmix_mag,
                    nelm=250,
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
                    setups={'base': 'recommended',
                            'W': '_sv'},
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
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True,
                    lsol=True,
                    nupdown=spin
                    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)