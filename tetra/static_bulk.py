import time
import subprocess
import numpy as np
from os import path
from mendeleev import element
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'static_bulk'
start_time = time.time()

spin_states_plus_2 = {'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 4,
                      'Co': 3, 'Ni': 2, 'Cu': 1, 'Zn': 1, 'Ga': 1, 'Ge': 2,
                      'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 4,
                      'Rh': 3, 'Pd': 2, 'Ag': 1, 'Cd': 1, 'In': 1, 'Sn': 2,
                      'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 4,
                      'Ir': 3, 'Pt': 2, 'Au': 1, 'Hg': 1, 'Tl': 1, 'Pb': 2
                      }

ldau_luj = {'Ti':{'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr':{'L':2, 'U':3.50,  'J':0.0},
            'Mn':{'L':2, 'U':3.75, 'J':0.0},
            'Fe':{'L':2, 'U':4.30,  'J':0.0},
            'Co':{'L':2, 'U':3.32, 'J':0.0},
            'Ni':{'L':2, 'U':6.45, 'J':0.0},
            'Cu':{'L':2, 'U':9.00,  'J':0.0}
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
else:
    raise ValueError('Missing: restart.json')

lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
        
def get_bands(atoms):
    nbands = 0
    nbands_per_orbital = {'s': 1, 'p': 3, 'd': 5, 'f': 7}
    for symbol in atoms.get_chemical_symbols():
        if symbol == 'H':  # H is bugged
            nbands += 1
            continue
        orbitals = element(symbol).ec.get_valence().to_str().split()
        nbands += sum(nbands_per_orbital[orbital[1]] * int(orbital[2:]) for orbital in orbitals)
    return nbands
    
def extract_kpoints(file_path):
    with open(file_path, 'r') as file:
        file_content = file.readlines()
    nkx, nky, nkz = map(int, file_content[3].split())
    return ((nkx, nky, nkz))

nbands = get_bands(atoms)
kpoints = extract_kpoints('./isif2/KPOINTS')

atoms.calc = vasp_calculator.Vasp(
                    encut=600,
                    xc='PBE',
                    gga='PE',
                    prec='Normal',
                    #inimix=0,
                    #amix=0.05,
                    #bmix=0.0001,
                    #amix_mag=0.05,
                    #bmix_mag=0.0001,
                    kpts=kpoints,
                    kpar=4,
                    npar=16,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    nelm=200,
                    algo='Normal',
                    isif=2,
                    nsw=200, ## lobster
                    ibrion=2,
                    ediff=1e-6,
                    ediffg=-0.02,
                    lreal='False',
                    lasph=True, 
                    lvtot=False,
                    laechg=True,
                    isym=0, ## lobster
                    ispin=2,
                    lorbit=11,
                    ldau=True,
                    ldautype=2,
                    ldau_luj=ldau_luj,
                    ldauprint=2,
                    lmaxmix=lmaxmix,
                    setups={'base': 'recommended', 'W': '_sv'},
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    # nupdown=0
                    )

energy = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("="*40 + "\n")

print(f"Execution time logged in 'time.log'.")