from ase.io import read, write
import sys

if len(sys.argv) != 2:
    print("Usage: python script.py <z_shift_value>")
    sys.exit(1)
    
dz = float(sys.argv[1])
index = 100
atoms = read('restart.json')
for atom in atoms:
    if atom.symbol not in ['C', 'N', 'O', 'H']:
        atom.z = 10.0 + 0.2 * dz
        index = atom.index
    elif atom.index == index+1 or atom.symbol == 'O':
        atom.z = 12.0 + 0.2 * dz
    elif atom.index == index+2 and atom.symbol == 'H':
        atom.z = 12.6 + 0.2 * dz
        
write('restart.json', atoms)