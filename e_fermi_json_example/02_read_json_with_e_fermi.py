#!/usr/bin/env python3
"""
Example 2: How to read JSON that contains e_fermi with ASE.

Reads sample.json produced by 01_add_e_fermi_to_json.py and
accesses e_fermi via atoms.info['e_fermi'].
"""
from read_json_with_info import read_json_with_info

# ---------- Configuration (JSON file to read) ----------
json_path = 'sample.json'
# ------------------------------------------------------

def main():
    atoms = read_json_with_info(json_path)
    e_fermi = atoms.info.get('e_fermi')
    print(f'File: {json_path}')
    print(f'atoms.info["e_fermi"] = {e_fermi}')
    print(f'Structure: {len(atoms)} atoms, formula = {atoms.get_chemical_formula()}')


if __name__ == '__main__':
    main()
