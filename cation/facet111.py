import sys
import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.geometry import get_duplicate_atoms
from ase.constraints import FixAtoms

# 1️⃣ Read input structure
input_file = sys.argv[1]
atoms = read(input_file)

# # 2️⃣ Remove duplicate atoms within 1 Å
# duplicates = get_duplicate_atoms(atoms, cutoff=1.0)
# atoms = atoms[[i for i in range(len(atoms)) if i not in duplicates]]

# 3️⃣ Remove the top 18 atoms based on the z-axis position
atoms = Atoms([atom for atom in sorted(atoms, key=lambda atom: atom.position[2], reverse=True)[18:]],
              cell=atoms.cell, pbc=atoms.pbc)

# 4️⃣ Shift all atoms by 0.2 Å in the z-direction
atoms.positions[:, 2] += 0.2

# 5️⃣ Wrap atoms inside the unit cell
atoms.wrap()

# 6️⃣ Fix the lowest two layers (18 atoms in total)
atoms_sorted = sorted(atoms, key=lambda atom: atom.position[2])  # Sort by z position
indices_to_fix = [atom.index for atom in atoms_sorted[:18]]  # First 18 atoms
atoms.set_constraint(FixAtoms(indices=indices_to_fix))

# 7️⃣ Compute slab height and add 18 Å vacuum
vacuum = 18.0
min_z = atoms.positions[:, 2].min()  # Lowest z position
max_z = atoms.positions[:, 2].max()  # Highest z position
height = max_z - min_z + vacuum  # Slab thickness + vacuum

# 8️⃣ Update the unit cell with new height
l1, l2 = atoms.cell.lengths()[:2]  # Keep original x and y lengths
a1, a2, a3 = atoms.cell.angles()  # Keep original angles
atoms.cell = (l1, l2, height, a1, a2, a3)

# 9️⃣ Generate output filename with .json extension
output_file = os.path.splitext(input_file)[0] + ".json"

# 🔟 Save the modified structure
write(output_file, atoms)

# 🔍 Print information for verification
print(f"Processed {input_file} → Saved as {output_file}")
print(f"Final number of atoms: {len(atoms)}")