import pandas as pd
import re

def extract_last_iteration_occupancies(outcar_path):
    with open(outcar_path, 'r') as file:
        lines = file.readlines()

    # Variables to store the data of the last iteration
    last_iteration_data = []
    is_last_iteration = False

    # Iterate through the lines and identify the start of iterations
    for line in lines:
        if "Iteration" in line:
            is_last_iteration = True
            last_iteration_data = []
        if is_last_iteration:
            last_iteration_data.append(line)

    # Extract occupancies for specified atoms (atom9 to atom16)
    atom_indices = range(9, 17)  # Atom indices from 9 to 16
    occupancies = {f"atom_{i}": [] for i in range(9, 17)}
    
    for atom_index in atom_indices:
        if atom_index < 10:
            atom_label = f"atom =   {atom_index}"
        else:
            atom_label = f"atom =  {atom_index}"
        for i, line in enumerate(last_iteration_data):
            if atom_label in line:
                for j in range(i, len(last_iteration_data)):
                    if "occupancies and eigenvectors" in last_iteration_data[j]:
                        occupancy_lines = last_iteration_data[j + 2: j + 12]
                        for k, occ_line in enumerate(occupancy_lines):
                            occupancy = float(occ_line.split()[2])
                            occupancies[f"atom_{atom_index}"].append(occupancy)
                        break
    return occupancies

# Extract occupancies
last_iteration_occupancies = extract_last_iteration_occupancies('OUTCAR')

# Create a DataFrame
df_last_iteration_occupancies = pd.DataFrame.from_dict(last_iteration_occupancies, orient='index', columns=[f"occ{i}" for i in range(1, 11)])

df_last_iteration_occupancies.to_csv('occ.tsv', sep='\t', float_format='%.4f')