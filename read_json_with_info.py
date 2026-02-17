#!/usr/bin/env python3
"""
ASE JSON DB 파일을 읽을 때 key_value_pairs를 atoms.info에 풀어서 넣어줌.
예: e_fermi를 DB에 넣어두었다면 read_json_with_info() 후 atoms.info['e_fermi']로 바로 접근.
"""
from ase.db import connect


def read_json_with_info(filename, index=1):
    """
    ASE JSON DB 파일에서 Atoms를 읽고, key_value_pairs를 atoms.info에 병합.

    ase.io.read()는 key_value를 넘겨주지 않으므로, e_fermi 등 커스텀 값을 쓰려면
    이 함수를 쓰거나 ase.db.connect + get_atoms(add_additional_information=True) 후
    atoms.info['key_value_pairs']['e_fermi']로 접근해야 함.

    사용 예:
        atoms = read_json_with_info('test.json')
        e_fermi = atoms.info.get('e_fermi')
    """
    with connect(filename, serial=True) as db:
        atoms = db.get_atoms(index, add_additional_information=True)
    kv = atoms.info.pop('key_value_pairs', None)
    if kv:
        atoms.info.update(kv)
    return atoms
