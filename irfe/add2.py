from ase.io import read, write
from ase import Atoms

original_atoms = read('o-top-layer.vasp')

atoms = original_atoms.copy()
gap = atoms[31].position - atoms[30].position
for i in range(8):
    atoms[-i-1].position += gap/2
write('o-brg-layer.vasp', atoms)

for i in range(8):
    atoms += Atoms('H', positions=[atoms[-8].position + (0.0, -0.8, 0.5)])
write('oh-brg-layer.vasp', atoms)

atoms = original_atoms.copy()
for i in [33, 34, 36, 39]:
    atoms += Atoms('H', positions=[atoms[i].position + (0.9, 0.0, 0.0)])
for i in [32, 35, 37, 38]:
    atoms += Atoms('H', positions=[atoms[i].position + (0.5, 0.8, 0.0)])
write('oh-top-layer.vasp', atoms)

original_atoms = read('CONTCAR')

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[(atoms[30].position + atoms[31].position)/2 + (0.0, 0.0, 1.8)])
write('o-brg-Ir.vasp', atoms)
atoms += Atoms('H', positions=[atoms[-8].position + (0.0, -0.8, 0.5)])
write('oh-brg-Ir.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[(atoms[28].position + atoms[29].position)/2 + (0.0, 0.0, 1.8)])
write('o-brg-M.vasp', atoms)
atoms += Atoms('H', positions=[atoms[-8].position + (0.0, -0.8, 0.5)])
write('oh-brg-M.vasp', atoms)
