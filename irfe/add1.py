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

        if ads == 'o' and site == 'top':
            atoms += Atoms('H', positions=[atoms[-1].position + (0.6, 0.8, 0.0)])
            write(f'{ads}h-{site}-M.vasp', atoms)
        elif ads == 'o' and site == 'hol':
            atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 0.0, 1.0)])
            write(f'{ads}h-{site}-M.vasp', atoms)

        # Layer adsorption
        atoms = original_atoms.copy()
        for i in ranges:
            atoms += Atoms(symbol, positions=[atoms[i].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-layer.vasp', atoms)

        if ads == 'o' and site == 'top':
            for i in range(8):
                atoms += Atoms('H', positions=[atoms[-8].position + (0.6, 0.8, 0.0)])
            write(f'{ads}h-{site}-layer.vasp', atoms)
        elif ads == 'o' and site == 'hol':
            for i in range(8):
                atoms += Atoms('H', positions=[atoms[-8].position + (0.0, 0.0, 1.0)])
            write(f'{ads}h-{site}-layer.vasp', atoms)