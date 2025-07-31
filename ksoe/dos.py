import csv
import numpy as np
import pandas as pd
import argparse
import os
import subprocess

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate d-band centers for multiple metal-alloy combinations.")
parser.add_argument("--output", type=str, default='all_d_band_centers', help="Output file name to save the results (e.g., 'all_d_band_centers.tsv').")
args = parser.parse_args()

# Define the directory structure
metal_types = ["2_Co", "3_Ni", "4_Cu"]
alloying_elements = ["01_Sc", "02_Ti", "03_V", "04_Cr", "05_Mn", "06_Fe", 
                     "07_Y", "08_Zr", "09_Nb", "10_Mo", "11_Ru", "12_W"]

# Define column names
columns = ["energy", "d(up)", "d(down)"]

# Function to create sumo.sh script for a specific metal
def create_sumo_script(directory, metal_name):
    sumo_script = f"""#!/bin/bash

sumo-dosplot \\
    --elements {metal_name}.d \\
    --atoms {metal_name}.19.20.22.23.25.26 \\
    --prefix surround

sumo-dosplot \\
    --elements {metal_name}.d \\
    --atoms {metal_name}.19.20.21.22.23.24.25.26 \\
    --prefix top

sumo-dosplot \\
    --elements {metal_name}.d \\
    --prefix host
"""
    
    script_path = os.path.join(directory, "sumo.sh")
    with open(script_path, 'w') as f:
        f.write(sumo_script)
    
    # Make the script executable
    os.chmod(script_path, 0o755)
    return script_path

# Function to calculate d-band center
def calculate_d_band_center(energy, dos_up, dos_down):
    dos_total = dos_up - dos_down
    numerator = np.trapz(energy * dos_total, energy)
    denominator = np.trapz(dos_total, energy)
    return numerator / denominator

# Store results
results = {}

# Process each combination
for metal_type in metal_types:
    # Extract metal name from metal_type (e.g., "2_Co" -> "Co")
    metal_name = metal_type.split('_')[1]
    
    # Define the DOS files for this metal type
    dos_files = [f"surround_{metal_name}_dos.dat", f"top_{metal_name}_dos.dat", f"host_{metal_name}_dos.dat"]
    
    for alloy_element in alloying_elements:
        directory = os.path.join(metal_type, alloy_element)
        
        # Initialize result entry for this combination
        key = (metal_type, alloy_element)
        results[key] = {"Surround": "N/A", "Top": "N/A", "Host": "N/A"}
        
        if os.path.exists(directory):
            # print(f"Processing {directory}...")
            
            # # Create and execute sumo script
            # try:
            #     script_path = create_sumo_script(directory, metal_name)
            #     print(f"  Running sumo.sh...")
                
            #     # Change to the directory and run the script
            #     result = subprocess.run(['bash', 'sumo.sh'], cwd=directory, 
            #                           capture_output=True, text=True, timeout=300)
                
            #     if result.returncode != 0:
            #         print(f"  Warning: sumo.sh execution had issues: {result.stderr}")
            #     else:
            #         print(f"  sumo.sh completed successfully")
                    
            # except subprocess.TimeoutExpired:
            #     print(f"  Error: sumo.sh execution timed out")
            # except Exception as e:
            #     print(f"  Error running sumo.sh: {e}")
            
            # Now process the DOS files
            for dos_file in dos_files:
                filepath = os.path.join(directory, dos_file)
                
                if os.path.exists(filepath):
                    try:
                        # Read the DOS data
                        data = pd.read_csv(filepath, sep='\s+', comment='#', names=columns)
                        
                        # Extract energy and DOS data
                        energy = data["energy"]
                        dos_up = data["d(up)"]
                        dos_down = data["d(down)"]
                        
                        # Calculate d-band center
                        d_band_center = calculate_d_band_center(energy, dos_up, dos_down)
                        
                        # Extract file type (surround, top, host)
                        file_type = dos_file.replace(f"_{metal_name}_dos.dat", "").replace("_", " ").title()
                        
                        results[key][file_type] = f"{d_band_center:.4f}"
                        print(f"  {file_type}: {d_band_center:.4f} eV")
                        
                    except Exception as e:
                        print(f"  Error processing {filepath}: {e}")
                        file_type = dos_file.replace(f"_{metal_name}_dos.dat", "").replace("_", " ").title()
                        results[key][file_type] = "Error"
                else:
                    print(f"  File not found: {filepath}")
                    file_type = dos_file.replace(f"_{metal_name}_dos.dat", "").replace("_", " ").title()
                    results[key][file_type] = "File not found"
        else:
            print(f"Directory not found: {directory}")
            results[key] = {"Surround": "Directory not found", "Top": "Directory not found", "Host": "Directory not found"}

# Save results to TSV file
tsv_file = f"{args.output}.tsv"
with open(tsv_file, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(["Host", "Dopant", "Surrnd", "Top", "Host"])
    for (metal_type, alloy_element), values in results.items():
        writer.writerow([metal_type, alloy_element, values["Surround"], values["Top"], values["Host"]])

print(f"\nResults saved as {tsv_file}")
print(f"Total processed: {len(results)} combinations")
