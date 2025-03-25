import pandas as pd
import subprocess
import os

# Define the full feature list
features = [
    'chg', 'mag', 'volume', 'l_bond', 'n_bond', 'grosspop', 'madelung', 
    'ICOHPm', 'ICOHPn', 'ICOBIm', 'ICOBIn', 'ICOOPm', 'ICOOPn', 
    'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
    'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 'Tboil', 'Tmelt', 
    'Hevap', 'Hfus', 'Hform'
]

# Initialize results DataFrames with rows as features and columns as iteration index
df_r2 = pd.DataFrame(index=features)
df_mae = pd.DataFrame(index=features)
df_mse = pd.DataFrame(index=features)

# Initial set of remaining features
remaining_features = features.copy()

# Iterate by progressively removing one feature at each step
for i in range(len(features)):
    for f in remaining_features:
        test_features = [x for x in remaining_features if x != f]
        input_str = ' '.join(test_features)
        output_name = f'seqleave_{i}_{f}'
        log_name = f'lr_{output_name}.log'

        # Run the regression script with the current feature set
        cmd = f'python ~/bin/tools/tetra/lr.py --Y form --X {input_str} --output {output_name}'
        subprocess.run(cmd, shell=True)

        # Parse the log file to extract performance metrics
        r2_val = mae_val = mse_val = None
        with open(log_name, 'r') as log_file:
            for line in log_file:
                if line.startswith('R2:'):
                    r2_val = float(line.split(':')[1].strip())
                elif line.startswith('MAE:'):
                    mae_val = float(line.split(':')[1].strip())
                elif line.startswith('MSE:'):
                    mse_val = float(line.split(':')[1].strip())

        df_r2.loc[f, i] = r2_val
        df_mae.loc[f, i] = mae_val
        df_mse.loc[f, i] = mse_val

    # Optional: choose a strategy for the next feature to remove
    # For now, just remove one (e.g., the one with highest MAE this round)
    next_feature = df_mae[i].astype(float).idxmax()  # or use df_r2[i].astype(float).idxmin()
    print(f"Iteration {i}: removing '{next_feature}'")
    remaining_features.remove(next_feature)

# Save results to files
df_r2.to_csv('recursive_leaveout_r2.tsv', sep='\t')
df_mae.to_csv('recursive_leaveout_mae.tsv', sep='\t')
df_mse.to_csv('recursive_leaveout_mse.tsv', sep='\t')
