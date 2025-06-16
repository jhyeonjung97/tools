import os
import csv

# Base directory for the oxidation state folders
base_dir = '/Users/jiuy97/Desktop/7_V_bulk/comer'
oxidation_dirs = {
    3: '1_Octahedral_+3',
    4: '2_Octahderal_+4',
    5: '3_Octahderal_+5',
    6: '4_Octahedral_+6'
}

def read_icohp_sum(filepath):
    """Read the -ICOHP sum from icohp-d.txt file"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            if lines:
                last_line = lines[-1].strip()
                if '-ICOHP sum:' in last_line:
                    return float(last_line.split(':')[1])
    except Exception as e:
        print(f"Error reading file {filepath}: {str(e)}")
        return None
    return None

# Read all_data.csv
csv_path = os.path.join(base_dir, 'all_data.csv')
print(f"Reading CSV file: {csv_path}")

with open(csv_path, 'r') as f:
    csv_reader = csv.DictReader(f)
    for row in csv_reader:
        try:
            metal = row['Metal Symbol']
            oxidation_state = int(float(row['Oxidation State']))
            bulk_icohp = abs(float(row['Bulk ICOHP']))
            
            print(f"\nProcessing {metal} (+{oxidation_state})")
            
            # Skip if oxidation state is not in our directories
            if oxidation_state not in oxidation_dirs:
                print(f"Skipping {metal} (+{oxidation_state}) - oxidation state not in directories")
                continue
            
            # Get metal number
            metal_num = None
            for i, m in enumerate(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'], 1):
                if m == metal:
                    metal_num = str(i).zfill(2)
                    break
            
            if metal_num is None:
                print(f"Skipping {metal} - metal not in list")
                continue
                
            # Construct path to icohp-d.txt
            metal_dir = f"{metal_num}_{metal}"
            icohp_path = os.path.join(base_dir, oxidation_dirs[oxidation_state], '3d', metal_dir, 'icohp-d.txt')
            print(f"Checking file: {icohp_path}")
            
            # Read icohp-d.txt value
            icohp_value = read_icohp_sum(icohp_path)
            
            if icohp_value is None:
                print(f"Could not read ICOHP value from {icohp_path}")
                continue
                
            print(f"CSV value: {bulk_icohp}, File value: {icohp_value}")
            
            # Compare values
            if abs(icohp_value - bulk_icohp) > 1e-5:
                print(f"Mismatch found for {metal} (+{oxidation_state}):")
                print(f"  all_data.csv value: {bulk_icohp}")
                print(f"  icohp-d.txt value: {icohp_value}")
                print(f"  Difference: {abs(icohp_value - bulk_icohp)}")
                print()
            else:
                print("Values match")
                
        except Exception as e:
            print(f"Error processing row: {row}")
            print(f"Error: {str(e)}")
            print() 