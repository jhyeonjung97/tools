import sys
import json
import subprocess
import numpy as np
from os import path
from math import sqrt
from mendeleev import element
from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'static_bulk2'

effective_length = 25

ldau_luj = {'Ti':{'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr':{'L':2, 'U':3.5,  'J':0.0},
            'Mn':{'L':2, 'U':3.75, 'J':0.0},
            'Fe':{'L':2, 'U':4.3,  'J':0.0},
            'Co':{'L':2, 'U':3.32, 'J':0.0},
            'Ni':{'L':2, 'U':6.45, 'J':0.0},
            'Cu':{'L':2, 'U':3.0,  'J':0.0},
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
else:
    raise ValueError('No restart.json file found')

lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
    
def get_bands(atoms):
    """
    Returns the exact number of bands desired by LOBSTER for the pCOHP calculations.
    """
    nbands = 0
    nbands_per_orbital = {'s': 1, 'p': 3, 'd': 5, 'f': 7}
    for symbol in atoms.get_chemical_symbols():
        if symbol == 'H':  # H is bugged
            nbands += 1
            continue
        orbitals = element(symbol).ec.get_valence().to_str().split()
        nbands += sum(nbands_per_orbital[orbital[1]] * int(orbital[2:]) for orbital in orbitals)
    return nbands

def get_kpoints(atoms, effective_length=effective_length, bulk=False):
    """
    Return a tuple of k-points derived from the unit cell.
    
    Parameters
    ----------
    atoms : object
    effective_length : k-point*unit-cell-parameter
    bulk : Whether it is a bulk system.
    """
    l = effective_length
    cell = atoms.get_cell()
    nkx = int(round(l/np.linalg.norm(cell[0]),0))
    nky = int(round(l/np.linalg.norm(cell[1]),0))
    if bulk == True:
        nkz = int(round(l/np.linalg.norm(cell[2]),0))
    else:
        nkz = 1
    return((nkx, nky, nkz))

nbands = get_bands(atoms)
kpoints = get_kpoints(atoms, effective_length=25, bulk=True)

atoms.calc = vasp_calculator.Vasp(
                    encut=600,
                    xc='PBE',
                    gga='PE',
                    kpts=kpoints,
                    kpar=8,
                    npar=1,
                    gamma=True,
                    ismear=0,
                    inimix=0,
                    amix=0.05,
                    bmix=0.0001,
                    amix_mag=0.01,
                    bmix_mag=0.00001,
                    nelm=250,
                    sigma=0.05,
                    algo='normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-6,
                    prec='Normal',
                    nsw=0,
                    lvtot=False,
                    nbands=nbands,
                    ispin=2,
                    setups={'base': 'recommended',
                            'W': '_sv'},
                    ldau=True,
                    ldautype=2,
                    laechg=True,
                    lreal='False',
                    lasph=True, 
                    ldau_luj=ldau_luj,
                    ldauprint=2,
                    lmaxmix=lmaxmix,
                    isym=0, 
                    nedos=3000,
                    lorbit=11
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    )

energy = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)