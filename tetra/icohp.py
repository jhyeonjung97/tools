"""
Module for reading ICOHPLIST.lobster and generating icohp.txt, icohp-s.txt, and icohp-d.txt
"""

import numpy as np
import pandas as pd
from ase.io import read
import re

def read_icohplist(filepath='.'):
    """
    Read ICOHPLIST.lobster file and return a pandas DataFrame
    """
    # Read ICOHPLIST.lobster
    with open(f'{filepath}/ICOHPLIST.lobster', 'r') as f:
        lines = f.readlines()
    
    # Read POSCAR to get element symbols
    atoms = read(f'{filepath}/POSCAR')
    
    # Initialize lists to store data
    data = []
    
    # Metal elements list
    metal_elements = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']
    
    # Dictionary to store ICOHP values for each bond
    bond_values = {}
    
    # Track current spin
    current_spin = None
    
    for line in lines:
        if line.strip() == '':
            continue
            
        # Split line into components
        parts = line.split()
        
        # Check for spin information
        if 'for spin' in line:
            current_spin = int(parts[-1])
            continue
            
        # Skip header line
        if 'COHP#' in line:
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
            
            # Get ICOHP value
            icohp = float(parts[7])
            
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
                    'total_icohp': 0.0,
                    'd_p_icohp': 0.0,
                    's_p_icohp': 0.0
                }
            
            # Update total ICOHP only if no orbitals are specified
            if orb1 is None and orb2 is None:
                bond_values[bond_id]['total_icohp'] += icohp
            
            # Update d-p ICOHP if metal d orbital and oxygen p orbital
            if (ele1 in metal_elements and orb1 is not None and 'd' in orb1 and 
                ele2 == 'O' and orb2 is not None and 'p' in orb2):
                bond_values[bond_id]['d_p_icohp'] += icohp
            elif (ele2 in metal_elements and orb2 is not None and 'd' in orb2 and 
                  ele1 == 'O' and orb1 is not None and 'p' in orb1):
                bond_values[bond_id]['d_p_icohp'] += icohp
            
            # Update s-p ICOHP if metal s orbital and oxygen p orbital
            if (ele1 in metal_elements and orb1 is not None and 's' in orb1 and 
                ele2 == 'O' and orb2 is not None and 'p' in orb2):
                bond_values[bond_id]['s_p_icohp'] += icohp
            elif (ele2 in metal_elements and orb2 is not None and 's' in orb2 and 
                  ele1 == 'O' and orb1 is not None and 'p' in orb1):
                bond_values[bond_id]['s_p_icohp'] += icohp
    
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
            'pair': f"{metal}{metal_idx}-{oxygen}{oxygen_idx}",
            'd_pair': f"{metal}{metal_idx}(d)-{oxygen}{oxygen_idx}(p)",
            's_pair': f"{metal}{metal_idx}(s)-{oxygen}{oxygen_idx}(p)",
            '-ICOHP': -values['total_icohp'],
            'd_p_ICOHP': -values['d_p_icohp'],
            's_p_ICOHP': -values['s_p_icohp'],
            'distance': values['distance']
        })
    
    if not data:
        raise ValueError("No ICOHP data found in the file. Please check the file format.")
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Sort by label
    df = df.sort_values('label')
    
    return df

def write_icohp_files(df, filepath='.'):
    """
    Write DataFrame to icohp.txt, icohp-d.txt, and icohp-s.txt files
    """
    # Format the DataFrame for output
    df['label'] = df['label'].astype(str).str.rjust(5)
    df['ele1'] = df['ele1'].astype(str).str.rjust(4)
    df['idx1'] = df['idx1'].astype(str).str.rjust(4)
    df['ele2'] = df['ele2'].astype(str).str.rjust(4)
    df['idx2'] = df['idx2'].astype(str).str.rjust(4)
    df['pair'] = df['pair'].astype(str).str.rjust(12)
    df['d_pair'] = df['d_pair'].astype(str).str.rjust(12)
    df['s_pair'] = df['s_pair'].astype(str).str.rjust(12)
    df['-ICOHP'] = df['-ICOHP'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['d_p_ICOHP'] = df['d_p_ICOHP'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['s_p_ICOHP'] = df['s_p_ICOHP'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['distance'] = df['distance'].astype(str)
    
    # Write header
    header = 'label ele1  idx1 ele2  idx2        pair  -ICOHP  distance\n'
    
    # Write total ICOHP file
    with open(f'{filepath}/icohp.txt', 'w') as f:
        f.write(header)
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['pair']} {row['-ICOHP']} {row['distance']}\n")
        f.write(f"\t\t\t-ICOHP sum:{df['-ICOHP'].astype(float).sum():.5f}\n")
    
    # Write d-p ICOHP file
    with open(f'{filepath}/icohp-d.txt', 'w') as f:
        f.write(header)
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['d_pair']} {row['d_p_ICOHP']} {row['distance']}\n")
        f.write(f"\t\t\t-ICOHP sum:{df['d_p_ICOHP'].astype(float).sum():.5f}\n")
    
    # Write s-p ICOHP file
    with open(f'{filepath}/icohp-s.txt', 'w') as f:
        f.write(header)
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['s_pair']} {row['s_p_ICOHP']} {row['distance']}\n")
        f.write(f"\t\t\t-ICOHP sum:{df['s_p_ICOHP'].astype(float).sum():.5f}\n")

def main():
    """
    Main function to read ICOHPLIST.lobster and write icohp files
    """
    try:
        df = read_icohplist()
        write_icohp_files(df)
        print("Successfully generated icohp.txt, icohp-d.txt, and icohp-s.txt")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 