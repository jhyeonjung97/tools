import os
import glob
import argparse
import numpy as np
from sys import argv
from ase.io import read, write, Trajectory

def main():

    parser = argparse.ArgumentParser(description='Process a water chain.')
    parser.add_argument('-n', type=int, default=8, help='Length of water chain')
    parser.add_argument('-z', action='store_false', default=True, help='Turn on "is_within_height"')
    parser.add_argument('file', type=str, help='Path to the input file')
    args = parser.parse_args()
    
    n = args.n
    z = args.z
    file = args.file
    atoms = read(file)

    name, ext = os.path.splitext(file)
    write(f'{name}_0{ext}', atoms)

    oxygen_indices = water_chain(atoms, n, z)

    hydrogen = None
    for i in range(len(oxygen_indices) - 1):
        oxygen1 = oxygen_indices[i]
        oxygen2 = oxygen_indices[i+1]
        hydrogen11 = None
        hydrogen12 = None
        hydrogen21 = None
        hydrogen22 = None
    
        for atom in atoms:
            if atom.symbol == 'H':
                d1 = np.linalg.norm(atoms[oxygen1].position - atom.position)
                d2 = np.linalg.norm(atoms[oxygen2].position - atom.position)
                if 0.5 < d1 < 1.5:
                    if atom.index == hydrogen:
                        continue
                    elif hydrogen11 is None:
                        hydrogen11 = atom.index
                    elif hydrogen12 is None:
                        hydrogen12 = atom.index
                    else:
                        print(f"Too many hydrogens {atom} were found near oxygen1: {oxygen1}")
                if 0.5 < d2 < 1.5:
                    if hydrogen21 is None:
                        hydrogen21 = atom.index
                    elif hydrogen22 is None:
                        hydrogen22 = atom.index
                    else:
                        print(f"Too many hydrogens {atom} were found near oxygen2: {oxygen2}")
        
        d1 = np.linalg.norm(atoms[hydrogen11].position - atoms[oxygen2].position)
        d2 = np.linalg.norm(atoms[hydrogen12].position - atoms[oxygen2].position)
        if d1 < d2:
            hydrogen = hydrogen11
        else:
            hydrogen = hydrogen12
    
        if hydrogen11 is None or hydrogen12 is None:
            print(f"Not all required hydrogens were found near the oxygen1: {oxygen1}")
            exit(1)
        if hydrogen21 is None or hydrogen22 is None:
            print(f"Not all required hydrogens were found near the oxygen2: {oxygen2}")
            exit(1)
                          
        O = atoms[oxygen2].position
        A1 = atoms[hydrogen21].position
        A2 = atoms[hydrogen22].position
        
        M = (A1 + A2) / 2
        B1 = np.cross(A1 - A2, O - M)
        B1 = B1 / np.linalg.norm(B1) * np.linalg.norm(A1 - M) + M
        B2 = 2 * M - B1
        C1 = 2 * O - B1
        C2 = 2 * O - B2
    
        d1 = np.linalg.norm(atoms[hydrogen].position - C1)
        d2 = np.linalg.norm(atoms[hydrogen].position - C2)
        if d1 < d2:
            atoms[hydrogen].position = C1
        else:
            atoms[hydrogen].position = C2
            
        write(f'{name}_{i+1}{ext}', atoms)
        
    check_with_carbons(name, ext, oxygen_indices)

def check_with_carbons(name, ext, oxygen_indices):
    combined_traj = Trajectory('water-chain.traj', 'w')
    for i in range(0, len(oxygen_indices)):
        file=f'{name}_{i}{ext}'
        new_name = os.path.splitext(file)[0] + '_with_carbon'
        atoms = read(file)
        for atom in atoms:
            if atom.index in oxygen_indices:
                atom.symbol = 'C'
        combined_traj.write(atoms)
    combined_traj.close()
            
def is_atom_in_cylinder(atom, z):
    cylinder_radius = 5  # Ångstrom
    cylinder_height = 16  # Ångstrom
    box_center_x = 15  # Ångstrom, half of 30 Å
    box_center_y = 15  # Ångstrom, half of 30 Å
    box_center_z = 20  # Ångstrom, half of 30 Å
    z_min = (box_center_z * 2 - cylinder_height) / 2  # Ångstrom, (40 - 20) / 2
    z_max = z_min + cylinder_height
    
    atom_position = atom.position
    
    distance_from_axis = np.linalg.norm(atom_position[:2] - np.array([box_center_x, box_center_y]))
    is_within_radius = distance_from_axis <= cylinder_radius
    is_within_height = z_min <= atom_position[2] <= z_max
    
    if z:
        return is_within_radius and is_within_height
    else:
        return is_within_radius

def water_chain(atoms, n, z):
    cation_found = False
    cations = ['Li', 'Na', 'K', 'Rb', 'Cs']
    for atom in atoms:
        if atom.symbol in cations:
            cation_found = True
            print(atom.symbol)
            oxygen_indices = [o.index for o in atoms if o.symbol == 'O' \
                              and 1.0 < np.linalg.norm(atom.position - o.position) < 3.0]

    if not cation_found:
        oxygen_indices = [o.index for o in atoms if o.symbol == 'O' and is_atom_in_cylinder(o, z)]

    chains = [[i] for i in oxygen_indices]  # Initialize chains with individual oxygen atoms
    new_chains = chains
    i = 0
    
    while new_chains and i < n:
        i += 1
        print(f'Search water chains of length {i}')
        chains = new_chains
        new_chains = []
        for chain in chains:
            last_oxygen_index = chain[-1]
            for h in atoms:
                if h.symbol == 'H' and is_atom_in_cylinder(h, z):
                    if 0.8 < np.linalg.norm(atoms[last_oxygen_index].position - h.position) < 1.2:
                        for o in atoms:
                            if o.symbol == 'O' and is_atom_in_cylinder(o, z) and o.index not in chain \
                            and 1.5 < np.linalg.norm(h.position - o.position) < 2.5:
                                new_chain = chain + [o.index]
                                new_chains.append(new_chain)

    d_max = 0
    for chain in chains:
        d = np.linalg.norm(atoms[chain[0]].position - atoms[chain[-1]].position)
        if d > d_max:
            d_max = d
            longest_chain = chain
            
    return longest_chain

if __name__ == '__main__':
    main()