#!/usr/bin/env python3
"""
Example 1: How to write e_fermi into a JSON file.

- Run this script as-is to create a test sample.json with e_fermi stored.
- For real use: set json_path to your file and set e_fermi from OUTCAR parsing, etc.
"""
from ase import Atoms
import ase.db

# ---------- Configuration (edit as needed) ----------
json_path = 'sample.json'
e_fermi = 5.123   # Example value; use value parsed from OUTCAR, etc. in practice
# ----------------------------------------------------

def make_sample_json_and_add_e_fermi():
    """Create a minimal test structure, write it to JSON, and store e_fermi."""
    # Simple structure (e.g. small RuO2-like cell)
    atoms = Atoms(
        'Ru2O4',
        positions=[[0, 0, 0], [1.5, 1.5, 1.5], [1.5, 1.5, 0], [1.5, 0, 1.5], [0, 1.5, 1.5], [0, 0, 1.5]],
        cell=[[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        pbc=True
    )
    with ase.db.connect(json_path, serial=True) as db:
        db.write(atoms, e_fermi=e_fermi)
    print(f'Saved: {json_path}  (e_fermi={e_fermi})')


def add_e_fermi_to_existing_json(path, e_fermi, row_id=1):
    """Add e_fermi to a specific row of an existing JSON file."""
    with ase.db.connect(path, serial=True) as db:
        db.update(row_id, e_fermi=e_fermi)
    print(f'Updated: {path}  row id={row_id}, e_fermi={e_fermi}')


if __name__ == '__main__':
    # Demo: create sample.json and store e_fermi
    make_sample_json_and_add_e_fermi()

    # To add e_fermi to an existing JSON only, use:
    add_e_fermi_to_existing_json('/Users/jiuy97/bin/tools/bulk-ruo2.json', e_fermi=1.9, row_id=1)
