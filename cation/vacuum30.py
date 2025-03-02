import os
import sys
from ase.io import read, write
import numpy as np

def process_structure(input_file):
    # Read the structure from VASP format
    atoms = read(input_file, format='vasp')
    waters = read('b0.vasp', format='vasp')

    if len(atoms) == 42:
        # Sort atoms by z-coordinate
        atoms_sorted = sorted(atoms, key=lambda atom: atom.position[2])
        
        # Remove the lowest 9 atoms
        del atoms[atoms_sorted[0].index:atoms_sorted[8].index+1]

        # Find the new lowest atom's z-coordinate
        z_min = min(atom.position[2] for atom in atoms)

        # Shift all atoms so that the lowest atom is at z = 0.2
        shift = 0.2 - z_min
        atoms.positions[:, 2] += shift

    # Find new zmin and zmax
    z_min_new = min(atom.position[2] for atom in atoms)
    z_max_new = max(atom.position[2] for atom in atoms)
    height = z_max_new - z_min_new + 35

    # Adjust vacuum to 30
    l1, l2 = atoms.cell.lengths()[:2]  # Keep original x and y lengths
    a1, a2, a3 = atoms.cell.angles()  # Keep original angles
    atoms.cell = (l1, l2, 40, a1, a2, a3)
    
    # water
    l1, l2 = waters.cell.lengths()[:2]  # Keep original x and y lengths
    a1, a2, a3 = waters.cell.angles()  # Keep original angles
    atoms.cell = (l1, l2, 40, a1, a2, a3)
    
    z_min = min(water.position[2] for water in waters)
    shift = z_max_new + 1.5 - z_min
    waters.positions[:, 2] += shift
    
    atoms += waters

    # Write the modified structure in VASP format
    output_file = os.path.splitext(input_file)[0] + ".traj"
    write(output_file, atoms)
    print(f"Modified structure saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py POSCAR")
        sys.exit(1)
    process_structure(sys.argv[1])
