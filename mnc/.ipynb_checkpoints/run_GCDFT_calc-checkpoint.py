#! /usr/bin/envpython 

### Created by: Md. Delowar Hossain ###

### Modified by: Evan Carlson ###
### Date: 08/25/2024 ###

### This script runs GCDFT calculations ###

### Imports ###
import os
import shutil
import subprocess
import numpy as np
from sys import path
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

## Set up the initial structure
name = 'solv_uncharged'
atoms=read('restart.json')
atoms.write(name+'_init'+'.traj')
atoms.write(name+'_init'+'.cif')

## Define K-points for the calculation
def get_kpoints(atoms, effective_length=25, bulk=False):
    """Return a tuple of k-points derived from the unit cell.

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

## Set the k-points
kpoints = get_kpoints(atoms, effective_length=20, bulk=False)


## Set up the VASP calculator in ASE
calc = vasp_calculator.Vasp(
            istart=1,
            encut=500,   # Planewave energy cutoff
            xc='PBE',    # Exchange-correlation functional
            gga='PE',    # GGA functional
            kpts = kpoints,   # k-points
            kpar = 5,    # k-point parallelization
            npar= 8,     # Parallelization
            gamma = True,     # Gamma-centered (defaults to Monkhorst-Pack)
            ismear=0,    # Methfessel-Paxton smearing
            inimix = 0,  # Initial mixing
            amix = 0.1,  # Mixing parameter
            bmix = 0.0001,    # Mixing parameter
            amix_mag= 0.1,    # Mixing parameter
            bmix_mag= 0.0001, # Mixing parameter
            nelm=250,    # Max number of electronic steps
            sigma = 0.05,     # Smearing width
            algo = 'normal',  # Electronic minimization algorithm
            ibrion=2,    # Ionic relaxation algorithm
            isif=2,      # Relaxation of ions and cell
            #nupdown=0,  # AFM vs FM
            ediffg=-0.02,     # Force threshold for ionic relaxation
            ediff=1e-5,  # Energy threshold for electronic minimization
            prec='normal',    # Precision
            nsw=500,     # Max number of ionic steps
            ispin=2,     # Spin-polarization 
            lsol=True,   # Solvation
            eb_k= 78.4,  # Solvent dielectric constant
            tau=0,       # Solvation cavitation
            lambda_d_k=3.000, # Solvent Debye screening length
            lvtot=False, 
            lvhar=True, 
            lwave=True,    # Write WAVECAR
            ldau=True,     # LDA+U
            ldautype=2,    # LDA+U type
            laechg=True,   # Write AECCAR0
            lasph =True,   # Write AECCAR2
            ldau_luj={
                      'Ti':{'L':2,  'U':3.00, 'J':0.0},
                      'Mn':{'L':2,  'U':2.75, 'J':0.0},  # Carlson
                      'O':{'L':-1,  'U':0.0,  'J':0.0}},
            ldauprint=2, 
            lmaxmix=6, 
            isym = 2,      # Symmetry
            lorbit=11,
            lreal  = 'Auto',  # Real-space projection
            #idipol=3,        # Dipole correction
            #dipol=(0.5, 0.5, 0.5),  # Dipole correction
            ldipol=False      # Dipole correction
            )

## Do the solvated, uncharged calculation
atoms.set_calculator(calc)
atoms.get_potential_energy()
atoms.get_forces()

## Write the final structure of the solvated, uncharged calculation
traj2=Trajectory('final_with_calculator.traj',  'w')  
traj2.write(atoms)

## Convert the final structure to JSON format
subprocess.call('ase convert -f final_with_calculator.traj  final_with_calculator.json', shell=True)
subprocess.call('ase convert -f final_with_calculator.json restart.json', shell=True)
subprocess.call('ase convert -f OUTCAR full_relax.json', shell=True)

## Do pDOS calculation on final uncharged structure
# subprocess.call('cp ~/bin/rapiDOS* .; python ./rapiDOS.py', shell=True)


## Determine the total number of electrons in the slab
total_electrons = calc.get_number_of_electrons() 


## Create a directory for each charge state
parent_dir = '/pscratch/sd/e/ezc/MnOOX_OER_project/surfaces/AlphaMnO2/110/Coverage5Layers/TopSites/OHBridge_OHTop/OTop/ExpSolvHex/HDownHOverO/VaspSol'
restart_path = parent_dir + '/restart.json'
wavecar_path = parent_dir + '/WAVECAR' 
run_vasp_path = parent_dir + '/run_vasp.py'
outcar_path = parent_dir + '/OUTCAR'


# Specify the charge states
fraction_charges = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0] ## Change this to the desired charge states

# Perform the calculations for each charge state
for fraction in fraction_charges:
   
   # Calculate the number of electrons for the charge state
   nelect = (total_electrons + fraction) 

   # Create a directory for the charge state and copy the necessary files
   folder_name = 'nelect_{}'.format(nelect)
   folder_path = os.path.join(parent_dir, folder_name)
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
   traj2=Trajectory('final_with_calculator.traj',  'w')
   traj2.write(new_atoms)
   subprocess.call('ase convert -f final_with_calculator.traj  final_with_calculator.json', shell=True)
   subprocess.call('ase convert -f final_with_calculator.json restart.json', shell=True)
   subprocess.call('ase convert -f OUTCAR full_relax.json', shell=True)
   subprocess.call('/global/cfs/cdirs/m2997/bin/get_restart4', shell=True)

   # Change back to the parent directory
   os.chdir('../')

# Create a directory for the uncharged calculation and copy the OUTCAR file
folder_name = 'nelect_{}'.format(total_electrons)
folder_path = os.path.join(parent_dir, folder_name)
os.makedirs(folder_path)
shutil.copy(outcar_path, folder_path)
