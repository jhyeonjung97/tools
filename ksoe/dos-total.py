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
metal_types = ["2_Co"]
alloying_elements = ["01_Sc", "02_Ti"]

# Define column names
columns = ["energy", "d(up)", "d(down)"]

# Function to create sumo.sh script for a specific metal
def create_sumo_script(directory, metal_name, alloy_element_name):
    sumo_script = f"""#!/bin/bash

sumo-dosplot \\
    --elements {alloy_element_name}.d \\
    --prefix dopant27-shift \\

sumo-dosplot \\
    --elements {alloy_element_name}.d \\
    --prefix dopant27-no-shift \\
    --no-shift

sumo-dosplot \\
    --elements {alloy_element_name}.d \\
    --prefix dopant27-zero-energy-0.1 \\
    --zero-energy 0.1

sumo-dosplot \\
    --elements {alloy_element_name}.d \\
    --prefix dopant27-zero-energy-0.2 \\
    --zero-energy 0.2
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
        alloy_element_name = alloy_element.split('_')[1]

        # Define the DOS file combinations with different shift options
        dos_combinations = {
            "brg_shift": [f"host19_{metal_name}_dos.dat", f"dopant27-shift_{alloy_element_name}_dos.dat"],
            "brg_no_shift": [f"host19_{metal_name}_dos.dat", f"dopant27-no-shift_{alloy_element_name}_dos.dat"],
            "brg_zero_01": [f"host19_{metal_name}_dos.dat", f"dopant27-zero-energy-0.1_{alloy_element_name}_dos.dat"],
            "brg_zero_02": [f"host19_{metal_name}_dos.dat", f"dopant27-zero-energy-0.2_{alloy_element_name}_dos.dat"],
        }

        # Initialize result entry for this combination
        key = (metal_type, alloy_element)
        
        # Create column names for all combinations and energy ranges
        all_combinations = ["brg_shift", "brg_no_shift", "brg_zero_01", "brg_zero_02"]
        energy_ranges = [
            ("full", None, None),           # -무한~+무한
            ("below0", None, 0.0),          # -무한~0
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
            
            # Create and execute sumo script
            try:
                script_path = create_sumo_script(directory, metal_name, alloy_element_name)
                print(f"  Running sumo.sh...")
                
                # Change to the directory and run the script
                result = subprocess.run(['bash', 'sumo.sh'], cwd=directory, 
                                      capture_output=True, text=True, timeout=300)
                
                if result.returncode != 0:
                    print(f"  Warning: sumo.sh execution had issues: {result.stderr}")
                else:
                    print(f"  sumo.sh completed successfully")
                    
            except subprocess.TimeoutExpired:
                print(f"  Error: sumo.sh execution timed out")
            except Exception as e:
                print(f"  Error running sumo.sh: {e}")
            
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
    all_combinations = ["brg_shift", "brg_no_shift", "brg_zero_01", "brg_zero_02"]
    energy_ranges = [
        ("full", None, None),           # -무한~+무한
        ("below0", None, 0.0),          # -무한~0
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
