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

name = 'static_bulk_afm_low'

effective_length = 25

spin_states_plus_4 = {'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 4,
                      'Co': 3, 'Ni': 2, 'Cu': 1, 'Zn': 1, 'Ga': 1, 'Ge': 2,
                      'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 4,
                      'Rh': 3, 'Pd': 2, 'Ag': 1, 'Cd': 1, 'In': 1, 'Sn': 2,
                      'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 4,
                      'Ir': 3, 'Pt': 2, 'Au': 1, 'Hg': 1, 'Tl': 1, 'Pb': 2,
                      }

ldau_luj = {'Ti':{'L':2,  'U':3.00, 'J':0.0},
            'V': {'L':2,  'U':3.25, 'J':0.0},
            'Cr':{'L':2,  'U':3.5,  'J':0.0},
            'Mn':{'L':2,  'U':3.75, 'J':0.0},
            'Fe':{'L':2,  'U':4.3,  'J':0.0},
            'Co':{'L':2,  'U':3.32, 'J':0.0},
            'Ni':{'L':2,  'U':6.45, 'J':0.0},
            'Cu':{'L':2, 'U':3.0,  'J':0.0},
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
    i = 0
    for a in atoms:
        if a.symbol in spin_states_plus_4:
            a.magmom = i #*spin_states_plus_4.get(a.symbol)
            i *= -1 # set AFM
else:
    raise ValueError('Where is restart.json')
    
for a in atoms:
    if a.symbol not in ldau_luj:
        ldau_luj[a.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}

def get_bands(atoms):
    """
    returns the extact number of bands desired by lobster for the pCOHP calculations
    """
    nbands = 0
    for sym in atoms.get_chemical_symbols():
        if sym == 'H': # H is bugged
            nbands += 1
            continue
        config = element(sym).ec.get_valence().to_str()
        config = config.split()
        for c in config:
            if 's' in c:
                nbands += 1*int(c[-1])
            elif 'p' in c:
                nbands += 3*int(c[-1])
            elif 'd' in c:
                nbands += 5*int(c[-1])
            elif 'f' in c:
                nbands += 7*int(c[-1])
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
                    # istart=0,
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
                    nelm=600,
                    sigma=0.05,
                    algo='normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-6,
                    prec='Normal',
                    nsw=0, # static
                    lvtot=False,
                    nbands=nbands, # static
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
                    isym=0, # static
                    nedos=3000,
                    lorbit=11,
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    nupdown=0
                    )

eng = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
# subprocess.call(f'ase convert -f final_{name}.traj restart.json', shell=True)
subprocess.call(f'cp OUTCAR OUTCAR_{name}', shell=True)