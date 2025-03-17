import os
import sys
import shutil
import subprocess
import numpy as np
from ase.io import read, Trajectory
from ase.calculators.vasp import Vasp

basename = os.path.basename(os.getcwd())

# Load atomic structure from available file
if os.path.exists('restart.json'):
    atoms = read('restart.json')
elif os.path.exists('start.traj'):
    atoms = read('start.traj')
else:
    raise ValueError("Neither restart.json nor start.traj file found")

# Define d-electron counts for transition metals
d_electrons = {
    'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 6, 'Co': 7, 'Ni': 8, 'Cu': 9, 'Zn': 10,
    'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 6, 'Rh': 7, 'Pd': 8, 'Ag': 9, 'Cd': 10,
    'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 6, 'Ir': 7, 'Pt': 8, 'Au': 9, 'Hg': 10,
}

# Determine oxidation state (assuming MN4C system)
count_o = len([atom for atom in atoms if atom.symbol == 'O'])
count_h = len([atom for atom in atoms if atom.symbol == 'H'])
oxi = count_o * 2 - count_h + 2

# Assign magnetic moments based on d-electron count
for atom in atoms:
    if atom.symbol in d_electrons:
        d_electron = d_electrons[atom.symbol] - oxi
        if basename.lower() == 'ls':  # Low spin state
            if d_electron < 0 or d_electron > 10:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in LS state.")
            spin = 1 if d_electron % 2 == 1 else 0
        elif basename.lower() == 'is':  # Intermediate spin state
            if d_electron not in {4, 5, 6}:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in IS state.")
            spin = 3 if d_electron == 5 else 2
        elif basename.lower() == 'hs':  # High spin state
            if d_electron not in {2, 3, 4, 5, 6, 7, 8}:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in HS state.")
            spin = d_electron if d_electron < 5 else 10 - d_electron
        else:
            sys.exit(f"Error: Invalid spin state. Choose from 'ls', 'is', or 'hs'.")
        atom.magmom = spin

# Set up VASP calculator
calc = Vasp(
    encut=500,             # Energy cutoff for plane wave basis
    xc='PBE',              # Exchange-correlation functional
    gga='PE',              # Use PBE functional
    prec='normal',         # Precision mode
    setups='recommended',  # Recommended PAW setups
    kpts=(5, 5, 1),        # k-point mesh
    gamma=True,            # Use Gamma-point only k-space integration
    kpar=5,                # Parallelization over k-points
    npar=4,                # Number of parallel processes
    # inimix = 0,            # Initial mixing
    # amix = 0.1,            # Mixing parameter
    # bmix = 0.0001,         # Mixing parameter
    # amix_mag= 0.1,         # Mixing parameter
    # bmix_mag= 0.0001,      # Mixing parameter
    algo='normal',         # Algorithm for electronic minimization
    nelm=250,              # Maximum number of electronic steps
    ediff=1e-5,            # Energy convergence criterion
    ismear=0,              # Smearing method (Gaussian smearing)
    sigma=0.05,            # Smearing width
    ibrion=2,              # Ionic relaxation method
    isif=2,                # Stress tensor and relaxation settings
    nsw=200,               # Maximum number of ionic steps
    ediffg=-0.02,          # Force convergence criterion
    ldau=True,             # Enable LDA+U correction
    ldautype=2,            # Type of LDA+U approach
    ldau_luj = {
        'Cr': {'L': 2, 'U': 3.50, 'J': 0.00},
        'Fe': {'L': 2, 'U': 4.30, 'J': 0.00},
        'Co': {'L': 2, 'U': 3.32, 'J': 0.00},
        'Mn': {'L': 2, 'U': 3.75, 'J': 0.00},
        'Ni': {'L': 2, 'U': 6.45, 'J': 0.00},
        'Ti': {'L': 2, 'U': 3.00, 'J': 0.00},
        'V':  {'L': 2, 'U': 3.25, 'J': 0.00},
        'Cu': {'L': 2, 'U': 3.00, 'J': 0.00},
    },                     # LDA+U parameters
    ldauprint=2,           # Print LDA+U information
    lmaxmix=4,             # Maximum l quantum number for charge mixing
    ispin=2,               # Spin-polarized calculation
    # nupdown=0,             # AFM vs FM
    lorbit=11,             # Projected density of states (DOS) output
    ldipol=True,           # Enable dipole correction
    idipol=3,              # Direction of dipole correction
    dipol=(0.5, 0.5, 0.5), # Dipole correction center
    lasph=True,            # Non-spherical contributions in PAW potentials
    laechg=True,           # Write AECCAR charge densities
    lreal='auto',          # Real-space projection (speed optimization)
    lsol=True,             # Enable implicit solvation model
    eb_k=78.4,             # Dielectric constant for implicit solvation model
    lambda_d_k=3.0,        # Debye screening length for implicit solvation
    tau=0,                 # Solvent relaxation time
    lvhar=True,            # Write Hartree potential
    lvtot=False,           # Do not write total potential
)

# Assign calculator to atoms
atoms.set_calculator(calc)

# Run calculation
total_energy = atoms.get_potential_energy()
forces = atoms.get_forces()

# Save results
Trajectory('final_with_calculator.traj', 'w').write(atoms)
subprocess.call('ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
subprocess.call('~/bin/get_restart3', shell=True)

if not os.path.exists('DONE'):
    sys.exit("Calculation not completed.")
    
total_electrons = calc.get_number_of_electrons() 

# Create a directory for each charge state
parent_dir = os.path.getcwd()
restart_path = os.path.join(parent_dir, 'restart.json')
wavecar_path = os.path.join(parent_dir, 'WAVECAR')

# Specify the charge states
fraction_charges = [-0.5, -1.0, -1.5, -2.0, -2.5, -3.0] ## Change this to the desired charge states

# Perform the calculations for each charge state
for fraction in fraction_charges:
   
    # Calculate the number of electrons for the charge state
    nelect = total_electrons + fraction
    
    # Create a directory for the charge state and copy the necessary files
    folder_name = f'nelect_{nelect}'
    folder_path = os.path.join(parent_dir, folder_name)
    os.makedirs(folder_path)
    shutil.copy2(restart_path, folder_path)  
    shutil.copy2(wavecar_path, folder_path)
    
    # Change to the directory for the charge state and set up the calculation
    os.chdir(folder_path) 
    new_atoms=read('restart.json')
    
    # Run the charged calculation
    calc.set(nelect=nelect)
    new_atoms.set_calculator(calc)
    new_atoms.get_potential_energy()
    new_atoms.get_forces()
    
    # Write the final structure of the calculation
    Trajectory('final_with_calculator.traj', 'w').write(new_atoms)
    subprocess.call('ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
    subprocess.call('~/bin/get_restart3', shell=True)

# Create a directory for the uncharged calculation and copy the OUTCAR file
folder_name = 'nelect_{}'.format(total_electrons)
folder_path = os.path.join(parent_dir, folder_name)
os.makedirs(folder_path)

# Copy all files from parent_dir to folder_path
for file_name in os.listdir(parent_dir):
    src_file = os.path.join(parent_dir, file_name)
    dest_file = os.path.join(folder_path, file_name)

    # Only copy files, skip directories
    if os.path.isfile(src_file):
        shutil.copy2(src_file, dest_file)