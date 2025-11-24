from ase.io import read

atoms = read('restart.json')
a = atoms.cell[0][0]*1.01
b = atoms.cell[1][1]*1.01
c = atoms.cell[2][2]
atoms.set_cell([a, b, c], scale_atoms=True)
atoms.write('restart.json')