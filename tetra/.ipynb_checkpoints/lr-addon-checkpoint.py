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

# Initialize results DataFrames
df_r2 = pd.DataFrame(index=features)
df_mae = pd.DataFrame(index=features)
df_mse = pd.DataFrame(index=features)

# Initialize selected feature list
selected_features = []

# Track unselected features
remaining_features = features.copy()

for i in range(len(features)):
    for f in remaining_features:
        test_features = selected_features + [f]
        input_str = ' '.join(test_features)

        # Remove previous log if exists
        if os.path.exists('lr_addon.log'):
            os.remove('lr_addon.log')

        # Run lr.py with selected features
        cmd = f'python ~/bin/tools/tetra/lr.py --Y form --X {input_str} --output addon'
        subprocess.run(cmd, shell=True, check=True)

        # Parse lr_addon.log
        r2_val = mae_val = mse_val = None
        if os.path.exists('lr_addon.log'):
            with open('lr_addon.log', 'r') as log_file:
                for line in log_file:
                    if line.startswith('R2:'):
                        r2_val = float(line.split(':')[1].strip())
                    elif line.startswith('MAE:'):
                        mae_val = float(line.split(':')[1].strip())
                    elif line.startswith('MSE:'):
                        mse_val = float(line.split(':')[1].strip())
        else:
            print(f"Warning: lr_addon.log not found after testing with only '{f}'")

        # Store results
        df_r2.loc[f, i] = r2_val
        df_mae.loc[f, i] = mae_val
        df_mse.loc[f, i] = mse_val

    # Choose the next best feature to add (e.g., lowest MAE)
    next_feature = df_mae[i].astype(float).idxmin()
    print(f"Iteration {i}: adding '{next_feature}'")
    selected_features.append(next_feature)
    remaining_features.remove(next_feature)

    # Save results
    df_r2.to_csv('addon_r2.tsv', sep='\t')
    df_mae.to_csv('addon_mae.tsv', sep='\t')
    df_mse.to_csv('addon_mse.tsv', sep='\t')
