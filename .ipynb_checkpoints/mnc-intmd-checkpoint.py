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

name = 'mnc-intmd'

effective_length = 25

spin_states_plus_4 = {'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 2, 'Mn': 3, 'Fe': 2,
                      'Co': 3, 'Ni': 2, 'Cu': 1, 'Zn': 1, 'Ga': 1, 'Ge': 2,
                      'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 2, 'Tc': 3, 'Ru': 2,
                      'Rh': 3, 'Pd': 2, 'Ag': 1, 'Cd': 1, 'In': 1, 'Sn': 2,
                      'La': 1, 'Hf': 2, 'Ta': 3, 'W': 2, 'Re': 3, 'Os': 2,
                      'Ir': 3, 'Pt': 2, 'Au': 1, 'Hg': 1, 'Tl': 1, 'Pb': 2,
                      }

ldau_luj = {'Ti': {'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr': {'L':2, 'U':3.5, 'J':0.0},
            'Mn': {'L':2, 'U':3.75, 'J':0.0},
            'Fe': {'L':2, 'U':4.3, 'J':0.0},
            'Co': {'L':2, 'U':3.32, 'J':0.0},
            'Ni': {'L':2, 'U':6.45, 'J':0.0},
            'Cu': {'L':2, 'U':3.0, 'J':0.0},
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
elif path.exists('start.traj'):
    atoms = read('start.traj')
else:
    raise ValueError('Where is start.traj')

for a in atoms:
    if a.symbol in spin_states_plus_4:
        a.magmom = i*spin_states_plus_4.get(a.symbol)
        
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

# nbands = get_bands(atoms)
# kpoints = get_kpoints(atoms, effective_length=25, bulk=False)

atoms.calc = vasp_calculator.Vasp(
                    istart=0,
                    encut=500,
                    xc='PBE',
                    gga='PE',
                    kpts=(5,5,1),
                    kpar=8,
                    npar=1,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    # inimix=0,
                    # amix=0.05,
                    # bmix=0.0001,
                    # amix_mag=0.05,
                    # bmix_mag=0.0001,
                    # nelm=600,
                    algo='normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-6,
                    prec='Normal',
                    nsw=200,
                    lvhar=True,
                    lvtot=False,
                    # nbands=nbands,
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
                    # isym=0, 
                    nedos=3000,
                    lorbit=11,
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    # nupdown=0
                    )

eng = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
# subprocess.call(f'ase convert -f final_{name}.traj restart.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)