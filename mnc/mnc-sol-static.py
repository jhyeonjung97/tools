import re
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

name = 'mnc-static'

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
else:
    raise ValueError('No restart.json file found')
        
lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    elif atom.symbol in ['C', 'N', 'O', 'H']:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
    else:
        ldau_luj[atom.symbol] = {'L': 2, 'U': 0.0, 'J': 0.0}   

def get_bands(atoms, zval_mapping):
    nbands = 0
    for symbol in atoms.get_chemical_symbols():
        nbands += zval_mapping(symbol)
    return nbands

def parse_zval_from_potcar(potcar_content):
    symbols = []
    zvals = []
    symbol_pattern = re.compile(r"TITEL\s*=\s*PAW_PBE\s*(\w+)")
    zval_pattern = re.compile(r"ZVAL\s*=\s*(\d+\.\d+)")

    lines = potcar_content.splitlines()
    for line in lines:
        match1 = symbol_pattern.search(line)
        if match1:
            symbols.append(match1.group(1))
        match2 = zval_pattern.search(line)
        if match2:
            zvals.append(int(float(match2.group(1))))

    if len(symbols) != len(zvals):
        raise ValueError("Mismatch between number of symbols and ZVALs in POTCAR content.")

    zval_mapping = dict(zip(symbols, zvals))
    return zval_mapping

with open("POTCAR", "r") as potcar_file:
    potcar_content = potcar_file.read()

zval_mapping = parse_zval_from_potcar(potcar_content)
nbands = get_bands(atoms, zval_mapping)

print(nbands)

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
                    #amix_mag=0.01,
                    #bmix_mag=0.00001,
                    nelm=600,
                    algo='Normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-5,
                    prec='Normal',
                    nsw=200,
                    lvhar=True,
                    lvtot=False,
                    nbands=nbands,
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
                    isym=0,
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