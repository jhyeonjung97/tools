from ase.io import read, write
from ase import Atoms

original_atoms = read('CONTCAR')

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[atoms[37].position + (0.0, 0.0, 1.4)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.5, 0.8, 0.0)])
write('o-ooh.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('H', positions=[atoms[37].position + (0.5, 0.8, 0.0)])
write('o-oh.vasp', atoms)

atoms = read('o-oh.vasp')
del atoms[34]
write('v-oh.vasp', atoms)

del original_atoms[34]
write('v-o.vasp', original_atoms)

del original_atoms[34]
write('v-v.vasp', original_atoms)
