import os
import glob
import re
from ase.io import read, write

# Define root directories and subdirectory lists
# root_dirs = glob.glob("/pscratch/sd/j/jiuy97/6_MNC/pourbaix/*_*/")
root_dirs = ["/pscratch/sd/j/jiuy97/6_MNC/pourbaix/1_Fe/",
             "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/2_Co/",
             "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/3_Mo/"]
main_dirs = ["clean", "h", "o", "oh", "oh-o", "oho", "oh-oh", "ohoh", "o-oh", "ooh"]
sub_dirs = ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"]

# Set the regular expression pattern to find energy values
e0_pattern = re.compile(r"E0=\s*(-?\.\d+E[+-]?\d+)")

def find_min_e_dir(main_dir, sub_dirs):
    min_e0 = None
    min_e0_dir = None 
    for sub_dir in sub_dirs:
        oszicar_path = os.path.join(main_dir, sub_dir, "OSZICAR")
        if os.path.isfile(oszicar_path):
            with open(oszicar_path, 'r') as file:
                for line in file:
                    match = e0_pattern.search(line)
                    if match:
                        e0_value = float(match.group(1))
                        if min_e0 is None or e0_value < min_e0:
                            min_e0 = e0_value
                            min_e0_dir = sub_dir
    return min_e0_dir

for root_dir in root_dirs:
    os.chdir(root_dir)
    numb, metal = os.path.basename(os.path.normpath(root_dir)).split('_', 1)
    
    for main_dir in main_dirs:
        adsorbate = os.path.basename(os.path.normpath(main_dir))
        most_stable_dir = 'most_stable'  # Name of the most stable directory

        # Define the path to the CONTcar file
        contcar_path = os.path.join(root_dir, main_dir, most_stable_dir, 'CONTCAR')
        
        if os.path.exists(contcar_path):
            atoms = read(contcar_path)
        else:
            # Search for the subdirectory with the minimum energy
            min_e0_dir = find_min_e_dir(os.path.join(root_dir, main_dir), sub_dirs)
            if min_e0_dir is not None:
                atoms = read(os.path.join(root_dir, main_dir, min_e0_dir, 'CONTCAR'))
            else:
                print(f"No valid CONTCAR file found in {main_dir}")
                continue  # Skip to the next iteration if invalid

        vasp_output_file = os.path.join(root_dir, main_dir, f'{numb}{metal}{adsorbate}.vasp')
        write(vasp_output_file, atoms)
        print(f"Saved as {vasp_output_file}")

        png_file = os.path.join(root_dir, main_dir, f'{numb}{metal}{adsorbate}.png')
        write(png_file, atoms, rotation='-90x, -90y, 0z', show_unit_cell=0)

        print(f"Saved as {png_file}")
