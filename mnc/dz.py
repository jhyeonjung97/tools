from ase.io import read, write
import argv

atoms = read('restart.json')
for atom in atoms:
    if atom.symbol not in ['C', 'N', 'O', 'H']:
        atom.z = 10.0 + 0.2 * argv[1]
write('restart.json', atoms)