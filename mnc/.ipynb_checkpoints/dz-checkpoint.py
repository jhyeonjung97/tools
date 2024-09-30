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
        zi = atom.z
        atom.z = 10.0 + 0.2 * dz
        zf = atom.z
        index = atom.index
        zdiff = zf - zi
        print(zdiff)
    if atom.index > index:
        atom.z += zdiff
        
write('restart.json', atoms)