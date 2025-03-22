import os
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

# Settings
vacuum = 18.0              # Vacuum thickness in the z-direction (Å)
fixed_atom_count = 8       # Number of atoms to fix at the bottom
# target_min_z = 0.2         # Target minimum z-position after shifting (Å)

# Process all .vasp files in the current directory
for filename in os.listdir("."):
    if filename.endswith(".vasp"):
        try:
            # 1️⃣ Read input structure
            atoms = read(filename)

            # 2️⃣ Get current z-range of atoms
            z_positions = atoms.positions[:, 2]
            min_z = z_positions.min()
            max_z = z_positions.max()
            slab_thickness = max_z - min_z

            # # 3️⃣ Shift the slab so that the bottom is at target_min_z
            # atoms.positions[:, 2] -= min_z             # Align bottom to z = 0
            # atoms.positions[:, 2] += target_min_z      # Shift bottom to z = 0.2 Å

            # 4️⃣ Set new cell height (slab thickness + vacuum)
            new_height = slab_thickness + vacuum
            cell = atoms.cell.array.copy()
            cell[2, :] = [0.0, 0.0, new_height]         # Modify z component of the cell
            atoms.set_cell(cell)
            atoms.center(axis=2)                        # Center the structure in z

            # 5️⃣ Fix the bottom N atoms (based on z-position)
            sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
            indices_to_fix = [atom.index for atom in sorted_atoms[:fixed_atom_count]]
            atoms.set_constraint(FixAtoms(indices=indices_to_fix))

            # 6️⃣ Save the modified structure as .json
            output_file = os.path.splitext(filename)[0] + ".traj"
            write(output_file, atoms)

            # 7️⃣ Print confirmation
            print(f"✅ {filename} → {output_file} | Slab height: {slab_thickness:.2f} Å | Fixed atoms: {indices_to_fix}")

        except Exception as e:
            print(f"❌ Failed to process {filename}: {e}")
