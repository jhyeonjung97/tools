from ase.io import read, write
import numpy as np
atoms = read('final_with_calculator.json')

l1 = atoms.cell.lengths()[0]
l2 = atoms.cell.lengths()[1]
displacement = [l1*(1/2)*1/2, l2*(1/2-4/9)*np.sqrt(3)/2, 0]
atoms.translate(displacement)
atoms.wrap()
write('centered.json', atoms)