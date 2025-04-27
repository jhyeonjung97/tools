from ase.io import read, write
from ase import Atoms

original_atoms = read('CONTCAR')

for ads in ['h', 'o']:
    symbol = ads.upper()

    for site in ['top', 'hol']:
        if site == 'top':
            range0 = 25
            ranges = range(24, 32)
            dz = 1.8
        elif site == 'hol':
            range0 = 9
            ranges = range(8, 16)
            dz = 5.8

        # Single atom adsorption
        atoms = original_atoms.copy()
        atoms += Atoms(symbol, positions=[atoms[range0].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-M.vasp', atoms)

        # Layer adsorption
        atoms = original_atoms.copy()
        for i in ranges:
            atoms += Atoms(symbol, positions=[atoms[i].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-layer.vasp', atoms)