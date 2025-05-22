import os
import glob
import argparse
import subprocess
import numpy as np
from ase.io import read
from statistics import mean

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, default='final_with_calculator.json')
parser.add_argument('-m', '--magnetic', action='store_true', default=False, help='Print magnetic moments')
parser.add_argument('-c', '--charge', action='store_true', default=False, help='Print Bader charges')
parser.add_argument('-p', '--cell', action='store_true', default=False, help='Print Cell Parameters')
parser.add_argument('-e', '--energy', action='store_true', default=False, help='Print total energy')
parser.add_argument('-a', '--atoms', action='store_true', default=False, help='Print chemical formula')
parser.add_argument('-f', '--force', action='store_true', default=False, help='Force to select all directories')
parser.add_argument('-b', '--beta', action='store_true', default=False)
parser.add_argument('-v', '--volume', action='store_true', default=False)
parser.add_argument('-hh', '--hydrogen', action='store_true', default=False)
parser.add_argument('-z', action='store_true', default=False)
parser.add_argument('-aa', action='store_true', default=False)
args = parser.parse_args()

class Colors:
    CYAN = '\033[34m'
    GREEN = '\033[38;5;22m'
    ORANGE = '\033[38;5;208m'
    MAGENTA = '\033[35m'
    RED = '\033[91m'
    RESET = '\033[0m'

if args.force:
    dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
else:
    dirs = [d for d in os.listdir('.') if os.path.isdir(d) and '_' in d]
    
sorted_dirs = sorted(dirs)
for dir in sorted_dirs:
    dir_path = os.path.join('.', dir)

    # pattern_A = os.path.join(dir_path, 'restart.json')

    # matching_files = []
    # for pattern in [pattern_A]: #, pattern_B, pattern_C]:
    #     matching_files.extend(glob.glob(pattern))
    #     if matching_files:
    #         break
    
    # atoms = None
    # for traj_file in matching_files:
    #     if os.path.exists(traj_file):
    #         atoms = read(traj_file)
    #         break
    
    path = os.path.join(dir_path, args.input)
    if os.path.exists(path):
        atoms = read(path)
        if args.magnetic:
            try:
                for atom in atoms:
                    if atom.symbol not in ['C', 'H', 'O', 'N']:
                        print(f"{Colors.CYAN}{dir}{Colors.RESET}", atom.symbol, atom.index, 
                              atoms.get_magnetic_moments()[atom.index])
            except Exception as e:
                print(f"{Colors.CYAN}{dir}{Colors.RESET}\t", 'Magnetic moments not available for this calculation.')
        if args.energy:
            print(f"{Colors.ORANGE}{dir}{Colors.RESET}\t", atoms.get_total_energy())
        if args.atoms:
            print(f"{Colors.MAGENTA}{dir}{Colors.RESET}\t", atoms.get_chemical_formula())
        if args.beta:
            print(f"{Colors.GREEN}{dir}{Colors.RESET}\t{atoms.cell.cellpar()[4]:.2f}")
        if args.cell:
            print(f"{Colors.GREEN}{dir}{Colors.RESET}\t", atoms.cell.cellpar())
        if args.hydrogen:
            atoms = read(os.path.join(dir_path, 'atoms_bader_charge.json'))
            charges = atoms.get_initial_charges()
            h_charge = mean([charges[i] for i in range(len(atoms)) if atoms[i].symbol == 'H'])
            print(f"{Colors.GREEN}{dir}{Colors.RESET}\t", h_charge)
        if args.z:
            zN = mean([atom.z for atom in atoms if atom.symbol == 'N'])
            zM = mean([atom.z for atom in atoms if atom.symbol != 'N' and atom.symbol != 'C' and atom.symbol != 'O' and atom.symbol != 'H'])
            dz = abs(zN - zM)
            print(f"{Colors.GREEN}{dir}{Colors.RESET}\t", dz)
            # print(f"{Colors.GREEN}{dir}{Colors.RESET}\t", atoms[0].z, atoms[1].z)
        if args.aa:
            print(f"{Colors.GREEN}{dir}{Colors.RESET}\t{atoms.cell.cellpar()[1]:.4f}\t{atoms.cell.cellpar()[2]:.4f}")