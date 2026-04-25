import argparse

from ase.io import read, write

parser = argparse.ArgumentParser()
parser.add_argument("--element", required=True, help="Substitution element symbol (e.g., Ru, Fe, Co)")
args = parser.parse_args()

target_element = args.element
atoms = read('restart.json')

candidate_indices = [i for i, atom in enumerate(atoms) if atom.symbol not in {target_element, 'O'}]
if not candidate_indices:
    raise ValueError(f"No atoms found with symbols other than {target_element} and O.")

target_idx = max(candidate_indices, key=lambda i: atoms.positions[i][2])

atoms[target_idx].symbol = target_element

magmoms = atoms.get_initial_magnetic_moments()
if target_idx == 45:
    source_idx = 44
elif target_idx == 19:
    source_idx = 16
else:
    raise ValueError(f"Unsupported target_idx for magmom mapping: {target_idx}")

magmoms[target_idx] = magmoms[source_idx]

atoms.set_initial_magnetic_moments(magmoms)

write('restart-m2ru.json', atoms)