#!/usr/bin/env python
import os
import glob
import numpy as np
import shutil
from ase.io import read

def compare_atomic_positions(base_path):
    """
    Compare atomic positions between CONTCAR and vib/POSCAR files in folders matching the specified path pattern.
    If positions are different, print the path.
    
    Args:
        base_path: base path pattern to search for
    """
    all_dirs = glob.glob(base_path)
    
    for dir_path in all_dirs:
        done_path = os.path.join(dir_path, 'DONE')
        contcar_path = os.path.join(dir_path, 'CONTCAR')
        poscar_path = os.path.join(dir_path, 'vib', 'POSCAR')
        
        if not os.path.exists(done_path):
            continue
        
        contcar_exists = os.path.exists(contcar_path)
        poscar_exists = os.path.exists(poscar_path)
        
        if not contcar_exists and not poscar_exists:
            print(f"both files not found: {dir_path}")
            continue
        elif not contcar_exists:
            print(f"CONTCAR file not found: {dir_path}")
            continue
        elif not poscar_exists:
            print(f"vib/POSCAR file not found: {dir_path}")
            continue
        
        try:
            contcar = read(contcar_path, format='vasp')
            poscar = read(poscar_path, format='vasp')
            
            contcar_positions = contcar.get_positions()
            poscar_positions = poscar.get_positions()
            
            if len(contcar_positions) != len(poscar_positions):
                print(f"number of atoms is different: {dir_path}")
                continue
            
            if not np.allclose(contcar_positions, poscar_positions, atol=1e-5):
                print(f"atomic positions are different: {dir_path}")
                vib_folder = os.path.join(dir_path, 'vib')
                if os.path.exists(vib_folder):
                    try:
                        shutil.rmtree(vib_folder)
                        print(f"deleted vib folder: {vib_folder}")
                    except Exception as e:
                        print(f"error deleting vib folder: {vib_folder} - {str(e)}")
        
        except Exception as e:
            print(f"Error: {dir_path} - {str(e)}")

if __name__ == "__main__":
    path_pattern = "/home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*"
    # path_pattern = "/home/hyeonjung/scratch/4_IrFe3/sgr/*_*_*/*_*"
    compare_atomic_positions(path_pattern)
