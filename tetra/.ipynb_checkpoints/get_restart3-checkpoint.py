#!/usr/bin/env python

import sys
import subprocess
import numpy as np
from ase.io import read, write
from ase.io.trajectory import Trajectory

def get_energy_from_outcar():
    """Extracts the energy value from OUTCAR."""
    energy = subprocess.check_output("grep py= OUTCAR | tail -1 | awk '{print $7}'", shell=True)
    return float(energy)

def process_trajectory(file):
    """Reads a trajectory file and returns the last atomic structure."""
    traj = Trajectory(file)
    return traj[-1]

def check_vasp_convergence():
    """Checks if the VASP output file contains the convergence message."""
    try:
        result = subprocess.check_output(
            "grep -F 'reached required accuracy' vasp.out", 
            shell=True
        )
        return bool(result.strip())  # Returns True if the message is found
    except subprocess.CalledProcessError:
        return False  # Returns False if the message is not found

def main():
    if len(sys.argv) > 1:
        file = sys.argv[1]
        atoms = process_trajectory(file)
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    else:
        file = 'OUTCAR'
        energy = get_energy_from_outcar()
        atoms = read(file)
        forces = atoms.get_forces()

    # Convert atomic structure to JSON using ASE
    subprocess.check_output(f"ase convert -f -n -1 {file} -o json moments.json", shell=True)

    # Save atomic structure to a trajectory file
    traj2 = Trajectory('moments.traj', 'w')
    traj2.write(atoms, energy=energy)

    # Analyze magnetic moments and forces
    moments = []
    moments2 = []
    moms = atoms.get_magnetic_moments()
    atoms.set_initial_magnetic_moments(moms)
    write('restart.json', atoms)

    total_force = 0.0
    max_force = 0.0

    if len(sys.argv) > 2:
        element = sys.argv[2]

    for a, atom in enumerate(atoms):
        if abs(moms[a]) > 0.1:
            moments.append([atom.symbol, a, moms[a]])

        if len(sys.argv) > 2 and atom.symbol == element:
            moments2.append(moms[a])

        force = np.linalg.norm(forces[a])
        total_force += force
        max_force = max(max_force, force)

    # Check if both convergence conditions are met
    vasp_converged = check_vasp_convergence()
    force_converged = max_force < 0.02

    if len(sys.argv) > 2:
        print(moments2)
    else:
        print(f'Moments from file {file}')
        for moment in moments:
            print(moment, ', ')

        print('Forces')
        print(f'Largest force: {max_force:.6f}, Total force: {total_force:.6f}')

        # Convergence is confirmed only if both conditions are satisfied
        if force_converged and vasp_converged:
            print('CONVERGED')
            subprocess.run(
                "echo $PWD ': CONVERGED with ENERGY %.6f and MAX FORCE of %.3f' > DONE" % (energy, max_force), 
                shell=True, check=True
            )


if __name__ == "__main__":
    main()
