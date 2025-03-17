import sys
import subprocess
import numpy as np
from os import path
from ase.io import read, Trajectory
from ase.calculators.vasp import Vasp

name = 'mnc-ls'

# Load atomic structure from available file
if path.exists('restart.json'):
    atoms = read('restart.json')
elif path.exists('start.traj'):
    atoms = read('start.traj')
else:
    raise ValueError('Neither restart.json nor start.traj file found')

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
        if 'ls' in name.lower():  # Low spin state
            if d_electron < 0 or d_electron > 10:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in LS state.")
            spin = 1 if d_electron % 2 == 1 else 0
        elif 'is' in name.lower():  # Intermediate spin state
            if d_electron not in {4, 5, 6}:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in IS state.")
            spin = 3 if d_electron == 5 else 2
        elif 'hs' in name.lower():  # High spin state
            if d_electron not in {2, 3, 4, 5, 6, 7, 8}:
                sys.exit(f"Error: Unexpected d_electron value ({d_electron}) for {atom.symbol} in HS state.")
            spin = d_electron if d_electron < 5 else 10 - d_electron
        else:
            sys.exit(f"Error: Invalid spin state '{spin_state}'. Choose from 'ls', 'is', or 'hs'.")
        atom.magmom = spin

# Define LDA+U parameters for specific elements
ldau_luj = {
    'Cr': {'L': 2, 'U': 3.50, 'J': 0.00},
    'Fe': {'L': 2, 'U': 4.30, 'J': 0.00},
    'Co': {'L': 2, 'U': 3.32, 'J': 0.00},
    'Mn': {'L': 2, 'U': 3.75, 'J': 0.00},
    'Ni': {'L': 2, 'U': 6.45, 'J': 0.00},
    'Ti': {'L': 2, 'U': 3.00, 'J': 0.00},
    'V':  {'L': 2, 'U': 3.25, 'J': 0.00},
    'Cu': {'L': 2, 'U': 3.00, 'J': 0.00},
}

# Set LDA+U parameters for other elements
lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    elif atom.symbol in ['C', 'N', 'O', 'H']:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}
    else:
        ldau_luj[atom.symbol] = {'L': 2, 'U': 0.0, 'J': 0.0}

# Set up VASP calculator
calc = Vasp(
    # Fundamental settings
    encut=500,  # Energy cutoff for plane wave basis
    xc='PBE',  # Exchange-correlation functional
    gga='PE',  # Use PBE functional
    prec='normal',  # Precision mode
    setups='recommended',  # Recommended PAW setups
    
    # K-point and parallelization
    kpts=(5, 5, 1),  # k-point mesh
    gamma=True,  # Use Gamma-point only k-space integration
    kpar=8,  # Parallelization over k-points
    npar=4,  # Number of parallel processes
    
    # Electronic convergence
    algo='normal',  # Algorithm for electronic minimization
    nelm=250,  # Maximum number of electronic steps
    ediff=1e-5,  # Energy convergence criterion
    ismear=0,  # Smearing method (Gaussian smearing)
    sigma=0.05,  # Smearing width
    
    # Ionic relaxation
    ibrion=2,  # Ionic relaxation method
    isif=2,  # Stress tensor and relaxation settings
    nsw=200,  # Maximum number of ionic steps
    ediffg=-0.02,  # Force convergence criterion
    
    # LDA+U settings
    ldau=True,  # Enable LDA+U correction
    ldautype=2,  # Type of LDA+U approach
    ldau_luj=ldau_luj,  # LDA+U parameters
    ldauprint=2,  # Print LDA+U information
    lmaxmix=lmaxmix,  # Maximum l quantum number for charge mixing
    
    # Spin and magnetism
    ispin=2,  # Spin-polarized calculation
    lorbit=11,  # Projected density of states (DOS) output
    
    # Dipole correction
    ldipol=True,  # Enable dipole correction
    idipol=3,  # Direction of dipole correction
    dipol=(0.5, 0.5, 0.5),  # Dipole correction center
    
    # Additional settings
    lasph=True,  # Non-spherical contributions in PAW potentials
    laechg=True,  # Write AECCAR charge densities
    lreal='auto',  # Real-space projection (speed optimization)
    
    # Implicit solvation model
    lsol=True,  # Enable implicit solvation model
    eb_k=78.4,  # Dielectric constant for implicit solvation model
    lambda_d_k=3.0,  # Debye screening length for implicit solvation
    tau=0,  # Solvent relaxation time
    
    # Output control
    lvhar=True,  # Write Hartree potential
    lvtot=False,  # Do not write total potential
)

# Assign calculator to atoms
atoms.set_calculator(calc)

# Run calculation
total_energy = atoms.get_potential_energy()
forces = atoms.get_forces()

# Save results
Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
subprocess.call(f'ase convert -f OUTCAR full_relax_{name}.json', shell=True)