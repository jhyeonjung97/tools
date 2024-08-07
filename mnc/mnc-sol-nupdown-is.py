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

name = 'mnc-nupdown-is'

effective_length = 25

spin_states_plus_2 = {'Cr': 2, 'Mn': 3, 'Fe': 2,
                      'Mo': 2, 'Tc': 3, 'Ru': 2,
                      'W': 2, 'Re': 3, 'Os': 2
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
    magmoms = atoms.get_magnetic_moments()
    atoms.set_initial_magnetic_moments(magmoms)
    for atom in atoms:
        if atom.symbol not in ['C', 'N', 'O', 'H']:
            if atom.symbol in spin_states_plus_2:
                spin = spin_states_plus_2.get(atom.symbol)
elif path.exists('start.traj'):
    atoms = read('start.traj')
    for atom in atoms:
        if atom.symbol not in ['C', 'N', 'O', 'H']:
            if atom.symbol in spin_states_plus_2:
                spin = spin_states_plus_2.get(atom.symbol)
                atom.magmom = spin
            else:
                raise ValueError(f"Unexpected atom symbol '{atom.symbol}' found in start.traj")
else:
    raise ValueError('Where is start.traj')
        
lmaxmix = 2
for atom in atoms:
    if atom.symbol in spin_states_plus_2:
        if atom.symbol in ldau_luj:
            lmaxmix = 4
        else:
            ldau_luj[atom.symbol] = {'L': 2, 'U': 0.0, 'J': 0.0}
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}

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
                    gga='PE',
                    ivdw=12,
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
                    algo='Normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-5,
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
                    # ldipol=True,
                    nupdown=spin,
                    lsol=True
                    )

e = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

traj_filename = f'final_{name}.traj'
Trajectory(traj_filename, 'w').write(atoms)
subprocess.call(f'ase convert -f {traj_filename} final_with_calculator.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)
