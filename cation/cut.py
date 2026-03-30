from ase.io import read, write

atoms = read('restart.json')

z_displacement = atoms.cell.lengths()[2] / 6 * 2
for atom in atoms:
    if atom.symbol == 'Li' or atom.symbol == 'Na':
        z_displacement = atoms.cell.lengths()[2] / 6 * 3        

atoms.positions[:, 2] -= z_displacement
atoms.wrap()

l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
l3 = atoms.cell.lengths()[2]
a1 = atoms.cell.angles()[0]
a2 = atoms.cell.angles()[1]
a3 = atoms.cell.angles()[2]
atoms.cell = (l1, l2, 40, a1, a2, a3)

# min_z = atoms.positions[:, 2].min()
# atoms.positions[:, 2] -= min_z - 10
atoms.center(axis=2)
write('restart_cut.json', atoms)