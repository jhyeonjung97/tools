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

name = 'opt_bulk3_n'

effective_length = 25

ldau_luj = {'Ti':{'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr':{'L':2, 'U':3.50, 'J':0.0},
            'Mn':{'L':2, 'U':3.75, 'J':0.0},
            'Fe':{'L':2, 'U':4.30, 'J':0.0},
            'Co':{'L':2, 'U':3.32, 'J':0.0},
            'Ni':{'L':2, 'U':2.00, 'J':0.0},
            'Cu':{'L':2, 'U':3.00, 'J':0.0},
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
elif path.exists('start.traj'):
    atoms = read('start.traj')
    for i in [8, 11, 13, 14]:
        atoms[i].magmom = 1
    for i in [9, 10, 12, 15]:
        atoms[i].magmom = -1
else:
    raise ValueError('Where is start.traj')
    
lmaxmix = 2
for a in atoms:
    if a.symbol in ldau_luj:
        lmaxmix = 4
    else:
        ldau_luj[a.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}

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

kpoints = get_kpoints(atoms, effective_length=effective_length, bulk=True)

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
                    amix_mag=0.05,
                    bmix_mag=0.0001,
                    nelm=250,
                    sigma=0.05,
                    algo='normal',
                    ibrion=2,
                    isif=3,
                    ediffg=-0.02,
                    ediff=1e-6,
                    prec='Normal',
                    nsw=600,
                    lvtot=False,
                    # nbands=nbands,
                    ispin=1,
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
                    # isym=0, 
                    nedos=3000,
                    lorbit=11,
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    nupdown=0
                    )

energy = atoms.get_potential_energy()
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)
