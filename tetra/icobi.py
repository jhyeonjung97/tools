"""
Module for reading ICOBILIST.lobster and generating icobi.txt and icobi-{d_orbital}2p.txt
"""

import numpy as np
import pandas as pd
from ase.io import read
import re
import os

def read_icobilist(filepath='.'):
    """
    Read ICOBILIST.lobster file and return a pandas DataFrame
    """
    # Read ICOBILIST.lobster
    with open(f'{filepath}/ICOBILIST.lobster', 'r') as f:
        lines = f.readlines()
    
    # Read atoms from final_with_calculator.json or moments.json
    if os.path.exists(f'{filepath}/final_with_calculator.json'):
        atoms = read(f'{filepath}/final_with_calculator.json')
    else:
        atoms = read(f'{filepath}/moments.json')
    
    # Initialize lists to store data
    data = []
    
    # Metal elements list with d orbital types
    metals = {
        '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
        '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
        '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
    }
    
    # Orbital order list (from lowest to highest energy)
    orbital_order = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p', '7s', '5f', '6d', '7p']
    
    # Dictionary to store ICOBI values for each bond
    bond_values = {}
    
    # Track current spin
    current_spin = None
    
    # Dictionary to store highest energy orbitals for each bond
    highest_orbitals = {}
    
    # Dictionary to store ICOBI values for each spin
    spin_values = {}
    
    # First pass: Find highest energy orbitals for each bond
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
        if 'COBI#' in line:
            continue
            
        # Check if this is a new bond
        if len(parts) >= 8 and parts[0].isdigit():
            bond_id = int(parts[0])
            atom1 = parts[1]
            atom2 = parts[2]
            
            # Get element symbols from atoms
            idx1 = int(re.findall(r'\d+', atom1)[0]) - 1  # Convert to 0-based index
            idx2 = int(re.findall(r'\d+', atom2)[0]) - 1
            ele1 = atoms[idx1].symbol
            ele2 = atoms[idx2].symbol
            
            # Get orbital information
            orb1 = parts[1].split('_')[1] if '_' in parts[1] else None
            orb2 = parts[2].split('_')[1] if '_' in parts[2] else None
            
            # Find highest energy orbitals
            if ele1 in metals['3d'] or ele1 in metals['4d'] or ele1 in metals['5d']:
                if orb1 is not None and orb2 is not None and '2p' in orb2:
                    if bond_id not in highest_orbitals:
                        highest_orbitals[bond_id] = {
                            'orb1': orb1,
                            'orb2': orb2,
                            'ele1': ele1,
                            'ele2': ele2,
                            'idx1': idx1,
                            'idx2': idx2
                        }
                    else:
                        current_orb = orb1.split('_')[0]
                        stored_orb = highest_orbitals[bond_id]['orb1'].split('_')[0]
                        if orbital_order.index(current_orb) > orbital_order.index(stored_orb):
                            highest_orbitals[bond_id] = {
                                'orb1': orb1,
                                'orb2': orb2,
                                'ele1': ele1,
                                'ele2': ele2,
                                'idx1': idx1,
                                'idx2': idx2
                            }
            elif ele2 in metals['3d'] or ele2 in metals['4d'] or ele2 in metals['5d']:
                if orb1 is not None and orb2 is not None and '2p' in orb1:
                    if bond_id not in highest_orbitals:
                        highest_orbitals[bond_id] = {
                            'orb1': orb2,
                            'orb2': orb1,
                            'ele1': ele2,
                            'ele2': ele1,
                            'idx1': idx2,
                            'idx2': idx1
                        }
                    else:
                        current_orb = orb2.split('_')[0]
                        stored_orb = highest_orbitals[bond_id]['orb1'].split('_')[0]
                        if orbital_order.index(current_orb) > orbital_order.index(stored_orb):
                            highest_orbitals[bond_id] = {
                                'orb1': orb2,
                                'orb2': orb1,
                                'ele1': ele2,
                                'ele2': ele1,
                                'idx1': idx2,
                                'idx2': idx1
                            }
    
    # Second pass: Read ICOBI values for the highest energy orbitals
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
        if 'COBI#' in line:
            continue
            
        # Check if this is a new bond
        if len(parts) >= 8 and parts[0].isdigit():
            bond_id = int(parts[0])
            atom1 = parts[1]
            atom2 = parts[2]
            distance = float(parts[3])
            
            # Get orbital information
            orb1 = parts[1].split('_')[1] if '_' in parts[1] else None
            orb2 = parts[2].split('_')[1] if '_' in parts[2] else None
            
            # Get ICOBI value
            icobi = float(parts[7])
            
            # Initialize bond_values if not exists
            if bond_id not in bond_values:
                if bond_id in highest_orbitals:
                    highest = highest_orbitals[bond_id]
                    bond_values[bond_id] = {
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'ele1': highest['ele1'],
                        'ele2': highest['ele2'],
                        'idx1': highest['idx1'],
                        'idx2': highest['idx2'],
                        'total_icobi': 0.0,
                        'd_p_icobi': 0.0,
                        'outer_icobi': 0.0,
                        'orb1': highest['orb1'],
                        'orb2': highest['orb2']
                    }
                    spin_values[bond_id] = {
                        'total_icobi': {1: 0.0, 2: 0.0},
                        'd_p_icobi': {1: 0.0, 2: 0.0},
                        'outer_icobi': {1: 0.0, 2: 0.0}
                    }
                else:
                    bond_values[bond_id] = {
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'ele1': atoms[int(re.findall(r'\d+', atom1)[0]) - 1].symbol,
                        'ele2': atoms[int(re.findall(r'\d+', atom2)[0]) - 1].symbol,
                        'idx1': int(re.findall(r'\d+', atom1)[0]) - 1,
                        'idx2': int(re.findall(r'\d+', atom2)[0]) - 1,
                        'total_icobi': 0.0,
                        'd_p_icobi': 0.0,
                        'outer_icobi': 0.0,
                        'orb1': None,
                        'orb2': None
                    }
                    spin_values[bond_id] = {
                        'total_icobi': {1: 0.0, 2: 0.0},
                        'd_p_icobi': {1: 0.0, 2: 0.0},
                        'outer_icobi': {1: 0.0, 2: 0.0}
                    }
            
            # Update total ICOBI only if no orbitals are specified
            if orb1 is None and orb2 is None:
                spin_values[bond_id]['total_icobi'][current_spin] += icobi
            
            # Update outer orbital ICOBI
            if bond_id in highest_orbitals:
                highest = highest_orbitals[bond_id]
                if orb1 == highest['orb1'] and orb2 == highest['orb2']:
                    spin_values[bond_id]['outer_icobi'][current_spin] += icobi
            
            # Update d-p ICOBI
            if bond_id in highest_orbitals:
                highest = highest_orbitals[bond_id]
                # Get the correct d orbital type for the metal
                metal_ele = highest['ele1']
                d_orbital = None
                if metal_ele in metals['3d']:
                    d_orbital = '3d'
                elif metal_ele in metals['4d']:
                    d_orbital = '4d'
                elif metal_ele in metals['5d']:
                    d_orbital = '5d'
                
                if d_orbital is not None:
                    if (str(orb1)==str(d_orbital) and str(orb2)=='2p') or (str(orb2)==str(d_orbital) and str(orb1)=='2p'):
                        spin_values[bond_id]['d_p_icobi'][current_spin] += icobi
    
    # Sum up spin values for each bond
    for bond_id in bond_values:
        if bond_id in spin_values:
            bond_values[bond_id]['total_icobi'] = sum(spin_values[bond_id]['total_icobi'].values())
            bond_values[bond_id]['d_p_icobi'] = sum(spin_values[bond_id]['d_p_icobi'].values())
            bond_values[bond_id]['outer_icobi'] = sum(spin_values[bond_id]['outer_icobi'].values())
    
    # Process bond_values to create data
    for bond_id, values in bond_values.items():
        # Determine which element is metal and which is oxygen
        if values['ele1'] in metals['3d'] or values['ele1'] in metals['4d'] or values['ele1'] in metals['5d']:
            metal, metal_idx = values['ele1'], values['idx1']
            oxygen, oxygen_idx = values['ele2'], values['idx2']
            orb_metal = values['orb1']
            orb_oxygen = values['orb2']
        else:
            metal, metal_idx = values['ele2'], values['idx2']
            oxygen, oxygen_idx = values['ele1'], values['idx1']
            orb_metal = values['orb2']
            orb_oxygen = values['orb1']
        
        data.append({
            'label': bond_id,
            'ele1': metal,
            'idx1': metal_idx,
            'ele2': oxygen,
            'idx2': oxygen_idx,
            'pair': f"{metal}{metal_idx}-{oxygen}{oxygen_idx}",
            'd_pair': f"{metal}{metal_idx}(d)-{oxygen}{oxygen_idx}(p)",
            'outer_pair': f"{metal}{metal_idx}({orb_metal})-{oxygen}{oxygen_idx}({orb_oxygen})",
            'ICOBI': values['total_icobi'],
            'd_p_ICOBI': values['d_p_icobi'],
            'outer_ICOBI': values['outer_icobi'],
            'distance': values['distance']
        })
    
    if not data:
        raise ValueError("No ICOBI data found in the file. Please check the file format.")
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Sort by label
    df = df.sort_values('label')
    
    return df

def write_icobi_files(df, filepath='.'):
    """
    Write DataFrame to icobi.txt, icobi-{d_orbital}2p.txt, and icobi-outer.txt files
    """
    # Metal elements list with d orbital types
    metals = {
        '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
        '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
        '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
    }
    
    # Format the DataFrame for output
    df['label'] = df['label'].astype(str).str.rjust(5)
    df['ele1'] = df['ele1'].astype(str).str.rjust(4)
    df['idx1'] = df['idx1'].astype(str).str.rjust(4)
    df['ele2'] = df['ele2'].astype(str).str.rjust(4)
    df['idx2'] = df['idx2'].astype(str).str.rjust(4)
    df['pair'] = df['pair'].astype(str).str.rjust(12)
    df['d_pair'] = df['d_pair'].astype(str).str.rjust(12)
    df['outer_pair'] = df['outer_pair'].astype(str).str.rjust(12)
    df['ICOBI'] = df['ICOBI'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['d_p_ICOBI'] = df['d_p_ICOBI'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['outer_ICOBI'] = df['outer_ICOBI'].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else x)
    df['distance'] = df['distance'].astype(str)
    
    # Write header
    header = 'label ele1  idx1 ele2  idx2        pair  ICOBI  distance\n'
    
    # Write total ICOBI file
    with open(f'{filepath}/icobi.txt', 'w') as f:
        f.write(header)
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['pair']} {row['ICOBI']} {row['distance']}\n")
        f.write(f"\t\t\tICOBI sum:{df['ICOBI'].astype(float).sum():.5f}\n")
    
    # Get d orbital type from the first metal element
    metal_ele = df['ele1'].iloc[0].strip()  # Remove whitespace
    d_orbital = None
    if metal_ele in metals['3d']:
        d_orbital = '3d'
    elif metal_ele in metals['4d']:
        d_orbital = '4d'
    elif metal_ele in metals['5d']:
        d_orbital = '5d'
    
    # Write d-p ICOBI file with dynamic filename
    if d_orbital is not None:
        with open(f'{filepath}/icobi-d.txt', 'w') as f:
            f.write(header)
            for _, row in df.iterrows():
                f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['d_pair']} {row['d_p_ICOBI']} {row['distance']}\n")
            f.write(f"\t\t\tICOBI sum:{df['d_p_ICOBI'].astype(float).sum():.5f}\n")
    
    # Write outer orbital ICOBI file
    with open(f'{filepath}/icobi-outer.txt', 'w') as f:
        f.write(header)
        for _, row in df.iterrows():
            f.write(f"{row['label']} {row['ele1']} {row['idx1']} {row['ele2']} {row['idx2']} {row['outer_pair']} {row['outer_ICOBI']} {row['distance']}\n")
        f.write(f"\t\t\tICOBI sum:{df['outer_ICOBI'].astype(float).sum():.5f}\n")
    
    return d_orbital

def main():
    """
    Main function to read ICOBILIST.lobster and write icobi files
    """
    try:
        # Check if ICOBILIST.lobster exists and unmatched file doesn't exist
        if not os.path.exists('ICOBILIST.lobster') or os.path.exists('unmatched'):
            return
            
        df = read_icobilist()
        d_name = write_icobi_files(df)
        print(f"Successfully generated icobi.txt, icobi-d.txt, and icobi-outer.txt")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 