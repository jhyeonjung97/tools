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

name = 'mnc-ls'

spin_states_plus_1 = {'Ti': 1, 'V': 0, 'Cr': 1, 'Mn': 0, 'Fe': 1, 'Co': 0, 'Ni': 1, 'Cu': 0,
                      'Zr': 1, 'Nb': 0, 'Mo': 1, 'Tc': 0, 'Ru': 1, 'Rh': 0, 'Pd': 1, 
                      'Hf': 1, 'Ta': 0, 'W': 1, 'Re': 0, 'Os': 1, 'Ir': 0, 'Pt': 1
                      }
spin_states_plus_2 = {'Ti': 0, 'V': 1, 'Cr': 0, 'Mn': 1, 'Fe': 0, 'Co': 1, 'Ni': 0, 'Cu': 1,
                      'Zr': 0, 'Nb': 1, 'Mo': 0, 'Tc': 1, 'Ru': 0, 'Rh': 1, 'Pd': 0, 
                      'Hf': 0, 'Ta': 1, 'W': 0, 'Re': 1, 'Os': 0, 'Ir': 1, 'Pt': 0
                      }
spin_states_plus_3 = {'Ti': 1, 'V': 0, 'Cr': 1, 'Mn': 0, 'Fe': 1, 'Co': 0, 'Ni': 1, 'Cu': 0,
                      'Zr': 1, 'Nb': 0, 'Mo': 1, 'Tc': 0, 'Ru': 1, 'Rh': 0, 'Pd': 1, 
                      'Hf': 1, 'Ta': 0, 'W': 1, 'Re': 0, 'Os': 1, 'Ir': 0, 'Pt': 1
                      }
spin_states_plus_4 = {'Ti': 0, 'V': 1, 'Cr': 0, 'Mn': 1, 'Fe': 0, 'Co': 1, 'Ni': 0, 'Cu': 1,
                      'Zr': 0, 'Nb': 1, 'Mo': 0, 'Tc': 1, 'Ru': 0, 'Rh': 1, 'Pd': 0, 
                      'Hf': 0, 'Ta': 1, 'W': 0, 'Re': 1, 'Os': 0, 'Ir': 1, 'Pt': 0
                      }
spin_states_plus_5 = {'Ti': 0, 'V': 0, 'Cr': 1, 'Mn': 0, 'Fe': 1, 'Co': 0, 'Ni': 1, 'Cu': 0,
                      'Zr': 0, 'Nb': 0, 'Mo': 1, 'Tc': 0, 'Ru': 1, 'Rh': 0, 'Pd': 1, 
                      'Hf': 0, 'Ta': 0, 'W': 1, 'Re': 0, 'Os': 1, 'Ir': 0, 'Pt': 1
                      }
spin_states_plus_6 = {'Ti': 0, 'V': 0, 'Cr': 0, 'Mn': 1, 'Fe': 0, 'Co': 1, 'Ni': 0, 'Cu': 1,
                      'Zr': 0, 'Nb': 0, 'Mo': 0, 'Tc': 1, 'Ru': 0, 'Rh': 1, 'Pd': 0, 
                      'Hf': 0, 'Ta': 0, 'W': 0, 'Re': 1, 'Os': 0, 'Ir': 1, 'Pt': 0
                      }
spin_states_plus_7 = {'Ti': 0, 'V': 0, 'Cr': 0, 'Mn': 0, 'Fe': 1, 'Co': 0, 'Ni': 1, 'Cu': 0,
                      'Zr': 0, 'Nb': 0, 'Mo': 0, 'Tc': 0, 'Ru': 1, 'Rh': 0, 'Pd': 1, 
                      'Hf': 0, 'Ta': 0, 'W': 0, 'Re': 0, 'Os': 1, 'Ir': 0, 'Pt': 1
                      }
spin_states_plus_8 = {'Ti': 0, 'V': 0, 'Cr': 0, 'Mn': 0, 'Fe': 0, 'Co': 1, 'Ni': 0, 'Cu': 1,
                      'Zr': 0, 'Nb': 0, 'Mo': 0, 'Tc': 0, 'Ru': 0, 'Rh': 1, 'Pd': 0, 
                      'Hf': 0, 'Ta': 0, 'W': 0, 'Re': 0, 'Os': 0, 'Ir': 1, 'Pt': 0
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
    #amix_mag = 0.01
    #bmix_mag = 0.00001
elif path.exists('start.traj'):
    atoms = read('start.traj')
    #amix_mag = 0.05
    #bmix_mag = 0.0001
else:
    raise ValueError('Neither restart.json nor start.traj file found')

count_c = len([atom for atom in atoms if atom.symbol == 'C']) - 26
count_n = len([atom for atom in atoms if atom.symbol == 'N'])
count_o = len([atom for atom in atoms if atom.symbol == 'O'])
count_h = len([atom for atom in atoms if atom.symbol == 'H'])
oxi = count_n / 2 + count_o * 2 - count_h

if count_c == 1 or oxi == 2:
    spin_states = spin_states_plus_2
elif oxi == 1:
    spin_states = spin_states_plus_1
elif oxi == 3:
    spin_states = spin_states_plus_3
elif oxi == 4:
    spin_states = spin_states_plus_4
elif oxi == 5:
    spin_states = spin_states_plus_5
elif oxi == 6:
    spin_states = spin_states_plus_6
elif oxi == 7:
    spin_states = spin_states_plus_7
elif oxi == 8:
    spin_states = spin_states_plus_8
else:
    spin_states = spin_states_plus_2
    
for atom in atoms:
    if atom.symbol in spin_states:
        spin = spin_states.get(atom.symbol)
        if not path.exists('restart.json'):
            atom.magmom = spin
    elif atom.symbol not in ['C', 'N', 'O', 'H']:
        raise ValueError(f"Unexpected atom symbol '{atom.symbol}' found in start.traj")
            
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
                    npar=4,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    #inimix=0,
                    #amix=0.05,
                    #bmix=0.0001,
                    #amix_mag=amix_mag,
                    #bmix_mag=bmix_mag,
                    nelm=250,
                    algo='Fast',
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
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True,
                    lsol=True
                    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)