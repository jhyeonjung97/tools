import os
import glob
from ase import Atoms
from ase.io import read, write

def main():
    # Loop through directories
    for dir in glob.glob('/pscratch/sd/j/jiuy97/4_V_slab/*_*_*/*d/*_*/'):
        # Skip directories with specific patterns in their names
        if 'x_' in dir or 's_' in dir or 'z_' in dir:
            continue
        if '8_Tetrahedral_AQ' in dir or '9_SquarePlanar_AU' in dir:
            continue
        
        # Change to the current directory
        os.chdir(dir)

        # Check if the 'DONE' file exists
        if not os.path.exists('DONE'):
            print(f"'DONE' file missing in directory: {dir}")
            continue

        # Read the atoms object from the file
        atoms = read('restart.json')
        
        # Assign marker and color based on the directory name and perform the respective operations
        metal_indices = [atom.index for atom in atoms if atom.symbol != 'O']
        
        if '1_Tetrahedral_WZ' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (-0.7, 0.0, 2.0)])
            write('restart-o.json', atoms)
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh.json', atoms)
            print(f'Saved json files in {dir}')
        
        elif '2_Tetrahedral_ZB' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-2]].position + (0.0, 0.0, 2.5)])
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o.json', atoms)
            atoms += Atoms('H', positions=[atoms[-2].position + (0.8, 0.0, 0.6)])
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh.json', atoms)
            print(f'Saved json files in {dir}')

        elif '3_SquarePlanar_TN' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-3]].position + (0.0, 0.0, 2.5)])
            atoms += Atoms('O', positions=[atoms[metal_indices[-2]].position + (0.0, 0.0, 2.5)])
            write('restart-o1.json', atoms)
            atoms += Atoms('H', positions=[atoms[-2].position + (0.8, 0.0, 0.6)])
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh1.json', atoms)
            # Reset atoms for the next operation
            atoms = read('restart.json')
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o2.json', atoms)
            atoms += Atoms('H', positions=[atoms[-2].position + (0.8, 0.0, 0.6)])
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh2.json', atoms)
            print(f'Saved json files in {dir}')

        elif '4_SquarePlanar_PD' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o1.json', atoms)
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh1.json', atoms)
            atoms = read('restart.json')
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o2.json', atoms)
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh2.json', atoms)
            print(f'Saved json files in {dir}')

        elif '5_SquarePlanar_NB' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-2]].position + (0.0, 0.0, 2.5)])
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o.json', atoms)
            atoms += Atoms('H', positions=[atoms[-2].position + (0.8, 0.0, 0.6)])
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh.json', atoms)
            print(f'Saved json files in {dir}')

        elif '6_Octahedral_RS' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-2]].position + (0.0, 0.0, 2.5)])
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o.json', atoms)
            atoms += Atoms('H', positions=[atoms[-2].position + (0.8, 0.0, 0.6)])
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh.json', atoms)
            print(f'Saved json files in {dir}')

        elif '7_Pyramidal_LT' in dir:
            atoms += Atoms('O', positions=[atoms[metal_indices[-1]].position + (0.0, 0.0, 2.5)])
            write('restart-o.json', atoms)
            atoms += Atoms('H', positions=[atoms[-1].position + (0.8, 0.0, 0.6)])
            write('restart-oh.json', atoms)
            print(f'Saved json files in {dir}')

if __name__ == '__main__':
    main()
