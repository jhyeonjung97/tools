import numpy as np
import pandas as pd
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate d-band centers for multiple metal-alloy combinations.")
parser.add_argument("--output", type=str, default='all_d_band_centers', help="Output file name to save the results (e.g., 'all_d_band_centers.tsv').")
args = parser.parse_args()

# Define column names
columns = ["energy", "d(up)", "d(down)"]

# Function to calculate d-band center
def calculate_d_band_center(energy, dos_up, dos_down):
    dos_total = dos_up - dos_down
    
    numerator = np.trapz(energy * dos_total, energy)
    denominator = np.trapz(dos_total, energy)
    return numerator / denominator

# Function to process DOS file and calculate band center
def process_dos_file(filename, element_name):
    try:
        if os.path.exists(filename):
            # Read the DOS data
            data = pd.read_csv(filename, sep='\s+', comment='#', names=columns)
            
            # Extract energy and DOS data
            energy = data["energy"]
            dos_up = data["d(up)"]
            dos_down = data["d(down)"]
            
            # Calculate d-band center
            d_band_center = calculate_d_band_center(energy, dos_up, dos_down)
            
            print(f"{element_name} band center: {d_band_center:.4f} eV")
            return d_band_center
        else:
            print(f"File not found: {filename}")
            return None
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

# Main execution
if __name__ == "__main__":
    # Process Ru DOS file
    ru_band_center = process_dos_file("sumo_Ru_dos.dat", "Ru")  
    # Process O DOS file
    o_band_center = process_dos_file("sumo_O_dos.dat", "O")
    
    # Save results to band-center.txt
    with open("band-center.txt", "w") as f:
        f.write("Band Center Results\n")
        f.write("==================\n\n")
        
        if ru_band_center is not None:
            f.write(f"Ru band center: {ru_band_center:.4f} eV\n")
        else:
            f.write("Ru band center: Error or file not found\n")
            
        if o_band_center is not None:
            f.write(f"O band center: {o_band_center:.4f} eV\n")
        else:
            f.write("O band center: Error or file not found\n")
    
    print(f"\nResults saved to band-center.txt")