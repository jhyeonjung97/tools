from ase.io import read, write
from ase.constraints import FixAtoms

atoms = read('restart.json')
l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
l3 = atoms.cell.lengths()[2]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.translate([l1/2, 0, 0])
atoms.cell = (l2, l3, l1, a2, a3, a1)
for atom in atoms:
    x = atom.position[0]
    y = atom.position[1]
    z = atom.position[2]
    atom.position = [y, z, x]
atoms = atoms.repeat((1, 1, 4))
sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
sorted_indices = [atom.index for atom in sorted_atoms]
re_index = min([atom.index for atom in sorted_atoms if atom.symbol != 'Ru'])
if atoms[re_index].symbol == 'O':
    re_index = 7
del_indices = []
for atom in atoms:
    if atom.position[2] < atoms[re_index].position[2]+1.5:
        del_indices.append(atom.index)
    elif atom.position[2] > atoms[re_index+72].position[2]-0.5:
        del_indices.append(atom.index)
for index in sorted(del_indices, reverse=True):
    del atoms[index]
if len(atoms) != 62:
    print('atoms is not 62')
    exit()
l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
l3 = atoms.cell.lengths()[2]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.cell = (l1, l2, 35, a1, a2, a3)
atoms.center(axis=2)
atoms.wrap()
sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
fixed = FixAtoms(indices=[atom.index for atom in sorted_atoms[:24]])
atoms.set_constraint(fixed)
write('o-sub-cus.json', atoms)

sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
del_indices = [sorted_atoms[-1].index, sorted_atoms[-2].index]
for index in sorted(del_indices, reverse=True):
    del atoms[index]
write('v-sub-cus.json', atoms)


atoms = read('restart.json')
l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
l3 = atoms.cell.lengths()[2]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.translate([l1/2, 0, 0])
atoms.cell = (l3, l1, l2, a3, a1, a2)
for atom in atoms:
    x = atom.position[0]
    y = atom.position[1]
    z = atom.position[2]
    atom.position = [z, x, y]
atoms = atoms.repeat((1, 1, 4))
sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
sorted_indices = [atom.index for atom in sorted_atoms]
re_index = min([atom.index for atom in sorted_atoms if atom.symbol != 'Ru'])
if atoms[re_index].symbol == 'O':
    re_index = 7
del_indices = []
for atom in atoms:
    if atom.position[2] < atoms[re_index].position[2]+1.5:
        del_indices.append(atom.index)
    elif atom.position[2] > atoms[re_index+72].position[2]-0.5:
        del_indices.append(atom.index)
for index in sorted(del_indices, reverse=True):
    del atoms[index]
if len(atoms) != 62:
    print('atoms is not 62')
    exit()
l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
l3 = atoms.cell.lengths()[2]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.cell = (l1, l2, 35, a1, a2, a3)
atoms.center(axis=2)
atoms.wrap()
sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
fixed = FixAtoms(indices=[atom.index for atom in sorted_atoms[:24]])
atoms.set_constraint(fixed)
write('o-sub-brg.json', atoms)

sorted_atoms = sorted(atoms, key=lambda atom: atom.position[2])
del_indices = [sorted_atoms[-1].index, sorted_atoms[-2].index]
for index in sorted(del_indices, reverse=True):
    del atoms[index]
write('v-sub-brg.json', atoms)
