import csv
import numpy as np
import pandas as pd
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate d-band centers for multiple metal-alloy combinations.")
parser.add_argument("--output", type=str, default='all_d_band_centers', help="Output file name to save the results (e.g., 'all_d_band_centers.tsv').")
args = parser.parse_args()

# Define the directory structure
metal_types = ["2_Co", "3_Ni", "4_Cu"]
alloying_elements = ["00_host", "01_Sc", "02_Ti", "03_V", "04_Cr", "05_Mn", "06_Fe", 
                     "07_Y", "08_Zr", "09_Nb", "10_Mo", "11_Ru", "12_W"]

# Define column names
columns = ["energy", "d(up)", "d(down)"]

# Function to create sumo.sh script for a specific metal
def create_sumo_script(directory, metal_name, alloy_element_name):
    sumo_script = f"""#!/bin/bash

sumo-dosplot \\
    --orbitals {metal_name}.d \\
    --elements {metal_name}.d \\
    --atoms {metal_name}.19 \\
    --prefix host19 \\

sumo-dosplot \\
    --orbitals {metal_name}.d \\
    --elements {metal_name}.d \\
    --atoms {metal_name}.20 \\
    --prefix host20 \\

sumo-dosplot \\
    --orbitals {alloy_element_name}.d \\
    --elements {alloy_element_name}.d \\
    --prefix dopant27 \\
"""
    
    script_path = os.path.join(directory, "sumo.sh")
    with open(script_path, 'w') as f:
        f.write(sumo_script)
    
    # Make the script executable
    os.chmod(script_path, 0o755)
    return script_path

# Function to calculate d-band center with optional energy range
def calculate_d_band_center(energy, dos_up, dos_down, energy_min=None, energy_max=None):
    dos_total = dos_up - dos_down
    
    # Apply energy filter if specified
    if energy_min is not None or energy_max is not None:
        mask = np.ones(len(energy), dtype=bool)
        if energy_min is not None:
            mask &= (energy >= energy_min)
        if energy_max is not None:
            mask &= (energy <= energy_max)
        
        energy_filtered = energy[mask]
        dos_total_filtered = dos_total[mask]
    else:
        energy_filtered = energy
        dos_total_filtered = dos_total
    
    numerator = np.trapz(energy_filtered * dos_total_filtered, energy_filtered)
    denominator = np.trapz(dos_total_filtered, energy_filtered)
    return numerator / denominator

# Store results
results = {}

# Process each combination
for metal_type in metal_types:
    # Extract metal name from metal_type (e.g., "2_Co" -> "Co")
    metal_name = metal_type.split('_')[1]
    
    for alloy_element in alloying_elements:
        directory = os.path.join(metal_type, alloy_element)
        if alloy_element == "00_host":
            alloy_element_name = metal_name
        else:
            alloy_element_name = alloy_element.split('_')[1]

        # Define the DOS file combinations
        dos_combinations = {
            "brg": [f"host19_{metal_name}_dos.dat", f"dopant27_{alloy_element_name}_dos.dat"],
            "hol": [f"host19_{metal_name}_dos.dat", f"host20_{metal_name}_dos.dat", f"dopant27_{alloy_element_name}_dos.dat"],
            "surround": [f"host19_{metal_name}_dos.dat", f"host20_{metal_name}_dos.dat", f"host22_{metal_name}_dos.dat", 
                       f"host23_{metal_name}_dos.dat", f"host25_{metal_name}_dos.dat", f"host26_{metal_name}_dos.dat", 
                       f"dopant27_{alloy_element_name}_dos.dat"],
            "layer": [f"host{i}_{metal_name}_dos.dat" for i in range(19, 27)] + [f"dopant27_{alloy_element_name}_dos.dat"]
        }

        # Initialize result entry for this combination
        key = (metal_type, alloy_element)
        
        # Create column names for all combinations and energy ranges
        all_combinations = ["brg", "hol", "surround", "layer"]
        energy_ranges = [
            ("full", None, None),           # -무한~+무한
            ("below0", None, 0.0),          # -무한~0
            ("m1to0", -1.0, 0.0),           # -1~0
            ("m2to0", -2.0, 0.0),           # -2~0
            ("m3to0", -3.0, 0.0),           # -3~0
            ("m4to0", -4.0, 0.0),           # -4~0
            ("m5to0", -5.0, 0.0),           # -5~0
            ("m6to0", -6.0, 0.0),           # -6~0
            ("m7to0", -7.0, 0.0),           # -7~0
            ("m8to0", -8.0, 0.0)            # -8~0
        ]
        
        results[key] = {}
        for comb in all_combinations:
            for range_name, _, _ in energy_ranges:
                if range_name == "full":
                    results[key][comb.title()] = "N/A"
                else:
                    results[key][f"{comb.title()}_{range_name}"] = "N/A"
        
        if os.path.exists(directory):
            print(f"Processing {directory}...")
            
            # Process each DOS combination
            for combination_name, dos_files in dos_combinations.items():
                try:
                    combined_energy = None
                    combined_dos_up = None
                    combined_dos_down = None
                    
                    files_found = 0
                    missing_files = []
                    
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
                                
                                if combined_energy is None:
                                    combined_energy = energy.copy()
                                    combined_dos_up = dos_up.copy()
                                    combined_dos_down = dos_down.copy()
                                else:
                                    # Add DOS data (assuming same energy grid)
                                    combined_dos_up += dos_up
                                    combined_dos_down += dos_down
                                
                                files_found += 1
                                
                            except Exception as e:
                                print(f"  Error reading {filepath}: {e}")
                                missing_files.append(dos_file)
                        else:
                            missing_files.append(dos_file)
                    
                    if files_found > 0 and combined_energy is not None:
                        # Calculate d-band center for all energy ranges
                        d_band_centers = []
                        for range_name, energy_min, energy_max in energy_ranges:
                            d_band_center = calculate_d_band_center(combined_energy, combined_dos_up, combined_dos_down, energy_min, energy_max)
                            d_band_centers.append(f"{d_band_center:.4f}")
                            
                            if range_name == "full":
                                results[key][combination_name.title()] = f"{d_band_center:.4f}"
                            else:
                                results[key][f"{combination_name.title()}_{range_name}"] = f"{d_band_center:.4f}"
                        
                        print(f"  {combination_name.title()}: {' | '.join(d_band_centers)} eV ({files_found}/{len(dos_files)} files)")
                        
                        if missing_files:
                            print(f"    Missing files: {', '.join(missing_files)}")
                    else:
                        for range_name, _, _ in energy_ranges:
                            if range_name == "full":
                                results[key][combination_name.title()] = "No valid files found"
                            else:
                                results[key][f"{combination_name.title()}_{range_name}"] = "No valid files found"
                        print(f"  {combination_name.title()}: No valid DOS files found")
                        
                except Exception as e:
                    print(f"  Error processing {combination_name}: {e}")
                    for range_name, _, _ in energy_ranges:
                        if range_name == "full":
                            results[key][combination_name.title()] = "Error"
                        else:
                            results[key][f"{combination_name.title()}_{range_name}"] = "Error"
                    
        else:
            print(f"Directory not found: {directory}")
            for comb in all_combinations:
                for range_name, _, _ in energy_ranges:
                    if range_name == "full":
                        results[key][comb.title()] = "Directory not found"
                    else:
                        results[key][f"{comb.title()}_{range_name}"] = "Directory not found"

# Save results to TSV file
tsv_file = f"{args.output}.tsv"
with open(tsv_file, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    
    # Create header with all combinations and energy ranges
    all_combinations = ["brg", "hol", "surround", "layer"]
    energy_ranges = [
        ("full", None, None),           # -무한~+무한
        ("below0", None, 0.0),          # -무한~0
        ("m1to0", -1.0, 0.0),           # -1~0
        ("m2to0", -2.0, 0.0),           # -2~0
        ("m3to0", -3.0, 0.0),           # -3~0
        ("m4to0", -4.0, 0.0),           # -4~0
        ("m5to0", -5.0, 0.0),           # -5~0
        ("m6to0", -6.0, 0.0),           # -6~0
        ("m7to0", -7.0, 0.0),           # -7~0
        ("m8to0", -8.0, 0.0)            # -8~0
    ]
    
    header = ["Host", "Dopant"]
    for comb in all_combinations:
        for range_name, _, _ in energy_ranges:
            if range_name == "full":
                header.append(comb.title())
            else:
                header.append(f"{comb.title()}_{range_name}")
    
    writer.writerow(header)
    
    for (metal_type, alloy_element), values in results.items():
        row = [metal_type, alloy_element]
        for comb in all_combinations:
            for range_name, _, _ in energy_ranges:
                if range_name == "full":
                    row.append(values[comb.title()])
                else:
                    row.append(values[f"{comb.title()}_{range_name}"])
        writer.writerow(row)

print(f"\nResults saved as {tsv_file}")
print(f"Total processed: {len(results)} combinations")
