"""
Module for reading ICOBILIST.lobster and generating icobi.txt
Hailey @SUNCAT/Stanford
Jul 1, 2023
"""

import numpy as np
import pandas as pd
from ase.io import read
import re

def read_icobilist(filepath='.'):
    """
    Read ICOBILIST.lobster file and return a pandas DataFrame
    """
    # Read ICOBILIST.lobster
    with open(f'{filepath}/ICOBILIST.lobster', 'r') as f:
        lines = f.readlines()
    
    # Read POSCAR to get element symbols
    atoms = read(f'{filepath}/POSCAR')
    
    # Initialize lists to store data
    data = []
    processed_bonds = set()  # To track processed bonds
    
    # Metal elements list from aloha/coop_analysis.py
    metal_elements = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']
    
    # Dictionary to store ICOBI values for each bond
    bond_values = {}
    
    for line in lines:
        if line.strip() == '':
            continue
            
        # Split line into components
        parts = line.split()
        
        # Skip header line
        if 'COBI#' in line:
            continue
            
        # Check if this is a new bond
        if len(parts) >= 8 and parts[0].isdigit():
            bond_id = int(parts[0])
            atom1 = parts[1]
            atom2 = parts[2]
            distance = float(parts[3])
            
            # Get element symbols from POSCAR
            idx1 = int(re.findall(r'\d+', atom1)[0]) - 1  # Convert to 0-based index
            idx2 = int(re.findall(r'\d+', atom2)[0]) - 1
            ele1 = atoms[idx1].symbol
            ele2 = atoms[idx2].symbol
            
            # Get orbital information
            orb1 = parts[1].split('_')[1] if '_' in parts[1] else None
            orb2 = parts[2].split('_')[1] if '_' in parts[2] else None
            
            if orb1 is None and orb2 is None:
                # This is a total ICOBI value
                icobi = float(parts[7])
                
                # Add to bond_values dictionary
                if bond_id not in bond_values:
                    bond_values[bond_id] = {
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'ele1': ele1,
                        'ele2': ele2,
                        'idx1': idx1,
                        'idx2': idx2,
                        'icobi': 0.0
                    }
                bond_values[bond_id]['icobi'] += icobi
    
    # Process bond_values to create data
    for bond_id, values in bond_values.items():
        # Determine which element is metal and which is oxygen
        if values['ele1'] in metal_elements:
            metal, metal_idx = values['ele1'], values['idx1']
            oxygen, oxygen_idx = values['ele2'], values['idx2']
        else:
            metal, metal_idx = values['ele2'], values['idx2']
            oxygen, oxygen_idx = values['ele1'], values['idx1']
        
        data.append({
            'label': bond_id,
            'ele1': metal,
            'idx1': metal_idx,
            'ele2': oxygen,
            'idx2': oxygen_idx,
            'pair': f"{metal}{metal_idx}(d)-{oxygen}{oxygen_idx}(p)",
            '-ICOBI': -values['icobi'],  # Change sign to match expected format
            'distance': values['distance']
        })
    
    if not data:
        raise ValueError("No ICOBI data found in the file. Please check the file format.")
    
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
        '-ICOBI': [df['-ICOBI'].sum()],
        'distance': ['']
    })
    df = pd.concat([df, sum_row], ignore_index=True)
    
    return df

def write_icobi_txt(df, filepath='.'):
    """
    Write DataFrame to icobi.txt file
    """
    # Format the DataFrame for output
    df['label'] = df['label'].astype(str).str.rjust(5)
    df['ele1'] = df['ele1'].astype(str).str.rjust(4)
    df['idx1'] = df['idx1'].astype(str).str.rjust(4)
    df['ele2'] = df['ele2'].astype(str).str.rjust(4)
    df['idx2'] = df['idx2'].astype(str).str.rjust(4)
    df['pair'] = df['pair'].astype(str).str.rjust(12)
    df['-ICOBI'] = df['-ICOBI'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['distance'] = df['distance'].astype(str)
    
    # Write header
    with open(f'{filepath}/icobi.txt', 'w') as f:
        f.write('label ele1  idx1 ele2  idx2         pair   -ICOBI  distance\n')
        
        # Write data
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['pair']} {row['-ICOBI']} {row['distance']}\n")

def main():
    """
    Main function to read ICOBILIST.lobster and write icobi.txt
    """
    try:
        df = read_icobilist()
        write_icobi_txt(df)
        print("Successfully generated icobi.txt")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 