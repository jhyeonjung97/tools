import sys
import json
import subprocess
import numpy as np
from os import path
from mendeleev import element
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'mnc-cep'

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

calc = vasp_calculator.Vasp(
    encut=500,
    xc='PBE',
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
    nelm=250,
    algo='Normal',
    prec='Normal',
    ibrion=2,
    isif=2,
    nsw=200,
    ediff=1e-5,
    ediffg=-0.02,
    setups={'base': 'recommended', 'W': '_sv'},
    lreal='Auto',
    lasph=True,
    laechg=True,
    ldau=True,
    ldautype=2,
    ldau_luj=ldau_luj,
    ldauprint=2,
    lmaxmix=lmaxmix,
    lorbit=11,
    ispin=2,
    ldipol=True,
    idipol=3,
    dipol=(0.5, 0.5, 0.5),
    lvhar=True,
    lvtot=False,
    lsol=True,
    eb_k=78.4,
    tau=0,
    lambda_d_k=3,
    nelect=146,
)

atoms.set_calculator(calc)
atoms.get_potential_energy()
atoms.get_forces()

print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)

total_electrons = calc.get_number_of_electrons() 

root = os.getcwd()
restart_path = root + '/restart.json'
wavecar_path = root + '/WAVECAR' 
run_vasp_path = root + '/run_vasp.py'
outcar_path = root + '/OUTCAR'


# Specify the charge states
fraction_charges = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0] ## Change this to the desired charge states

# Perform the calculations for each charge state
for fraction in fraction_charges:
   
   # Calculate the number of electrons for the charge state
   nelect = (total_electrons + fraction) 

   # Create a directory for the charge state and copy the necessary files
   folder_name = 'nelect_{}'.format(nelect)
   folder_path = os.path.join(root, folder_name)
   os.makedirs(folder_path)
   shutil.copy(restart_path, folder_path)  
   shutil.copy(wavecar_path, folder_path)
   shutil.copy(run_vasp_path, folder_path)

   # Change to the directory for the charge state and set up the calculation
   os.chdir(folder_path) 
   name = 'test_relax'
   new_atoms=read('restart.json')
   new_atoms.write(name+'_init'+'.traj')
   new_atoms.write(name+'_init'+'.cif')

   # Run the charged calculation
   calc.set(istart=1, nelect=nelect)
   new_atoms.set_calculator(calc)
   new_atoms.get_potential_energy()
   new_atoms.get_forces()

   # Write the final structure of the calculation
   traj2 = Trajectory('final_with_calculator.traj',  'w')
   traj2.write(new_atoms)
   subprocess.call('ase convert -f final_with_calculator.traj  final_with_calculator.json', shell=True)
   subprocess.call('ase convert -f final_with_calculator.json restart.json', shell=True)
   subprocess.call('ase convert -f OUTCAR full_relax.json', shell=True)
   subprocess.call('/global/cfs/cdirs/m2997/bin/get_restart4', shell=True)

   # Change back to the parent directory
   os.chdir('../')

# Create a directory for the uncharged calculation and copy the OUTCAR file
folder_name = 'nelect_{}'.format(total_electrons)
folder_path = os.path.join(root, folder_name)
os.makedirs(folder_path)
shutil.copy(outcar_path, folder_path)