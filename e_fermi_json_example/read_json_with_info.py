#!/usr/bin/env python3
"""
Read ASE JSON DB files and merge key_value_pairs into atoms.info.
If e_fermi (or other keys) are stored in the DB, use read_json_with_info()
and then access them as atoms.info['e_fermi'].
"""
from ase.db import connect


def read_json_with_info(filename, index=1):
    """
    Read Atoms from an ASE JSON DB file and merge key_value_pairs into atoms.info.

    Returns:
        Atoms: with e_fermi, etc. from key_value_pairs available in atoms.info.
    """
    with connect(filename, serial=True) as db:
        atoms = db.get_atoms(index, add_additional_information=True)
    kv = atoms.info.pop('key_value_pairs', None)
    if kv:
        atoms.info.update(kv)
    return atoms
