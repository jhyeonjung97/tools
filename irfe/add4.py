from ase.io import read, write
from ase import Atoms

original_atoms = read('CONTCAR')

for ads in ['h', 'o']:
    symbol = ads.upper()

    for site in ['top', 'hol']:
        if site == 'top':
            range1 = 7
            range2 = 11
            ranges = [6, 7, 8, 9, 10, 11, 12, 13]
            dz = 1.8
        elif site == 'hol':
            range1 = 2
            range2 = 23
            ranges = [2, 3, 20, 21,22, 23, 24, 25]
            dz = 5.8

        # Single atom adsorption
        atoms = original_atoms.copy()
        atoms += Atoms(symbol, positions=[atoms[range1].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-atom1.vasp', atoms)

        # Single atom adsorption
        atoms = original_atoms.copy()
        atoms += Atoms(symbol, positions=[atoms[range2].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-atom2.vasp', atoms)
        
        # Layer adsorption
        atoms = original_atoms.copy()
        for i in ranges:
            atoms += Atoms(symbol, positions=[atoms[i].position + (0.0, 0.0, dz)])
        write(f'{ads}-{site}-layer.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[(atoms[8].position + atoms[9].position)/2 + (0.0, 0.0, 1.8)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.0, -0.8, 0.5)])
write('oh-brg-atom1.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[(atoms[12].position + atoms[13].position)/2 + (0.0, 0.0, 1.8)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.0, -0.8, 0.5)])
write('oh-brg-atom2.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[atoms[7].position + (0.0, 0.0, 1.8)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 0.9, 0.0)])
write('oh-top-atom1.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[atoms[11].position + (0.0, 0.0, 1.8)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.0, 0.9, 0.0)])
write('oh-top-atom2.vasp', atoms)

original_atoms = read('o-top-layer.vasp')

atoms = original_atoms.copy()
gap = atoms[6].position - atoms[7].position
for i in range(8):
    atoms[-i-1].position += gap/2
for i in range(8):
    atoms += Atoms('H', positions=[atoms[-8].position + (0.0, -0.8, 0.5)])
write('oh-brg-layer.vasp', atoms)

atoms = original_atoms.copy()
for i in [33, 34, 36, 39]:
    atoms += Atoms('H', positions=[atoms[i].position + (0.9, 0.0, 0.0)])
for i in [32, 35, 37, 38]:
    atoms += Atoms('H', positions=[atoms[i].position + (0.5, 0.8, 0.0)])
write('oh-top-layer.vasp', atoms)
