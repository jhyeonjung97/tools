from ase.io import read, write
atoms=read('final_with_calculator.json')
for atom in atoms:
    if atom.symbol not in ['N', 'C', 'O', 'H']:
        dz = atom.z - 10.0
        mag = magmoms[atom.index]
        if atom.index < len(atoms) - 1:
            next_atom = atoms[atom.index + 1]
            print(atom.symbol, atom.position)
            print(next_atom.symbol, next_atom.position)
            bond = np.linalg.norm(atom.position - next_atom.position)
            print(bond)