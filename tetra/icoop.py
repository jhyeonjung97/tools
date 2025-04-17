"""
Module for reading ICOOPLIST.lobster and generating icoop.txt
Hailey @SUNCAT/Stanford
Jul 1, 2023
"""

import numpy as np
import pandas as pd
from ase.io import read
import re

def read_icooplist(filepath='.'):
    """
    Read ICOOPLIST.lobster file and return a pandas DataFrame
    """
    # Read ICOOPLIST.lobster
    with open(f'{filepath}/ICOOPLIST.lobster', 'r') as f:
        lines = f.readlines()
    
    # Read POSCAR to get element symbols
    atoms = read(f'{filepath}/POSCAR')
    
    # Initialize lists to store data
    data = []
    processed_bonds = set()  # To track processed bonds
    
    for line in lines:
        if line.strip() == '':
            continue
            
        # Split line into components
        parts = line.split()
        
        # Skip header line
        if 'COOP#' in line:
            continue
            
        # Check if this is a new bond
        if len(parts) >= 8 and parts[0].isdigit():
            bond_id = int(parts[0])
            atom1 = parts[1]
            atom2 = parts[2]
            distance = float(parts[3])
            
            # Skip if we've already processed this bond
            if bond_id in processed_bonds:
                continue
            
            # Get element symbols from POSCAR
            idx1 = int(re.findall(r'\d+', atom1)[0]) - 1  # Convert to 0-based index
            idx2 = int(re.findall(r'\d+', atom2)[0]) - 1
            ele1 = atoms[idx1].symbol
            ele2 = atoms[idx2].symbol
            
            # Get orbital information
            orb1 = parts[1].split('_')[1] if '_' in parts[1] else None
            orb2 = parts[2].split('_')[1] if '_' in parts[2] else None
            
            if orb1 is None and orb2 is None:
                # This is a total ICOOP value
                icoop = float(parts[7])
                
                # Determine which element is metal and which is oxygen
                if ele1 in ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']:
                    metal, metal_idx = ele1, idx1
                    oxygen, oxygen_idx = ele2, idx2
                else:
                    metal, metal_idx = ele2, idx2
                    oxygen, oxygen_idx = ele1, idx1
                
                # Always put metal as ele1
                data.append({
                    'label': bond_id,
                    'ele1': metal,
                    'idx1': metal_idx,
                    'ele2': oxygen,
                    'idx2': oxygen_idx,
                    'pair': f"{metal}{metal_idx}(d)-{oxygen}{oxygen_idx}(p)",
                    '-ICOOP': -icoop,  # Change sign to match expected format
                    'distance': distance
                })
                processed_bonds.add(bond_id)
    
    if not data:
        raise ValueError("No ICOOP data found in the file. Please check the file format.")
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Sort by label
    df = df.sort_values('label')
    
    # Add sum row
    sum_row = pd.DataFrame({
        'label': [''],
        'ele1': [''],
        'idx1': [''],
        'ele2': [''],
        'idx2': [''],
        'pair': [''],
        '-ICOOP': [df['-ICOOP'].sum()],
        'distance': ['']
    })
    df = pd.concat([df, sum_row], ignore_index=True)
    
    return df

def write_icoop_txt(df, filepath='.'):
    """
    Write DataFrame to icoop.txt file
    """
    # Format the DataFrame for output
    df['label'] = df['label'].astype(str).str.rjust(5)
    df['ele1'] = df['ele1'].astype(str).str.rjust(4)
    df['idx1'] = df['idx1'].astype(str).str.rjust(4)
    df['ele2'] = df['ele2'].astype(str).str.rjust(4)
    df['idx2'] = df['idx2'].astype(str).str.rjust(4)
    df['pair'] = df['pair'].astype(str).str.rjust(12)
    df['-ICOOP'] = df['-ICOOP'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['distance'] = df['distance'].astype(str)
    
    # Write header
    with open(f'{filepath}/icoop.txt', 'w') as f:
        f.write('label ele1  idx1 ele2  idx2         pair   -ICOOP  distance\n')
        
        # Write data
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['pair']} {row['-ICOOP']} {row['distance']}\n")

def main():
    """
    Main function to read ICOOPLIST.lobster and write icoop.txt
    """
    try:
        df = read_icooplist()
        write_icoop_txt(df)
        print("Successfully generated icoop.txt")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 