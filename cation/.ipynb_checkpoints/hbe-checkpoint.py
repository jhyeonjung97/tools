import os
from ase.io import read, write
from ase import Atoms, Atom

def add_hydrogen_and_save(atoms, indices, height, output_filename):
    """Add H atoms above specified indices at a given height and save the modified structure."""
    modified_atoms = atoms.copy()  # Copy original structure
    for index in indices:
        pos = modified_atoms[index].position.copy()  # Get original atom position
        pos[2] += height  # Move H atom in the z-direction
        modified_atoms.append(Atom("H", position=pos))  # ✅ Correct way to add H atom

    # Save to output JSON
    write(output_filename, modified_atoms)
    print(f"Saved {output_filename} with {len(indices)} added H atoms.")


# List of input files
input_files = ["a1.json", "a2.json", "a3.json", "a4.json"]

# Define target atom indices and heights for each modification
targets = [
    ([8], 1.8, "b"),  # bX.json
    ([15], 3.6, "c"),  # cX.json
    ([20], 6.0, "d"),  # dX.json
    ([0, 2, 4, 5, 6, 7], 1.8, "e"),  # eX.json
    ([9, 10, 16, 17, 14, 12], 3.6, "f"),  # fX.json
    ([23, 21, 19, 26, 25, 18], 6.0, "g"),  # hX.json
    (list(range(0, 9)), 1.8, "h"),  # iX.json (0~8)
    (list(range(9, 18)), 3.6, "i"),  # jX.json (9~17)
    (list(range(18, 27)), 6.0, "j")  # gX.json (18~26)
]

# Process each input file
for input_file in input_files:
    atoms = read(input_file)  # Read input structure
    prefix = os.path.splitext(input_file)[0][1]  # Extract number from "aX.json"

    # Apply modifications and save each version
    for indices, height, label in targets:
        output_file = f"{label}{prefix}.json"  # Generate output filename (e.g., b1.json, c1.json, ..., g4.json)
        add_hydrogen_and_save(atoms, indices, height, output_file)
