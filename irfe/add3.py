from ase.io import read, write
from ase import Atoms

original_atoms = read('o-o.vasp')

atoms = original_atoms.copy()
atoms += Atoms('O', positions=[atoms[33].position + (0.0, 0.0, 1.4)])
atoms += Atoms('H', positions=[atoms[-1].position + (0.5, 0.8, 0.0)])
write('o-ooh.vasp', atoms)

atoms = original_atoms.copy()
atoms += Atoms('H', positions=[atoms[33].position + (0.5, 0.8, 0.0)])
write('o-oh.vasp', atoms)

del original_atoms[35]
write('v-o.vasp', original_atoms)

del original_atoms[33]
write('v-v.vasp', original_atoms)

atoms = read('o-oh.vasp')
del atoms[35]
write('v-oh.vasp', atoms)

# mv CONTCAR 4_O_O
# mv o-oh.vasp 5_O_OH
# mv o-ooh.vasp 6_O_OOH
# mv v-oh.vasp 3_V_OH
# mv v-o.vasp 2_V_O
# mv v-v.vasp 1_V_V
