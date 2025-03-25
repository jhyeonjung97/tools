import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

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

# Initialize list of remaining features
remaining_features = features.copy()

for i in range(len(features)):
    for f in remaining_features:
        test_features = [x for x in remaining_features if x != f]
        input_str = ' '.join(test_features)

        # Remove previous log if exists
        if os.path.exists('lr_leaveout.log'):
            os.remove('lr_leaveout.log')

        # Run lr.py with selected features
        cmd = f'python ~/bin/tools/tetra/lr.py --Y form --X {input_str} --output leaveout'
        subprocess.run(cmd, shell=True, check=True)

        # Parse lr_leaveout.log
        r2_val = mae_val = mse_val = None
        if os.path.exists('lr_leaveout.log'):
            with open('lr_leaveout.log', 'r') as log_file:
                for line in log_file:
                    if line.startswith('R2:'):
                        r2_val = float(line.split(':')[1].strip())
                    elif line.startswith('MAE:'):
                        mae_val = float(line.split(':')[1].strip())
                    elif line.startswith('MSE:'):
                        mse_val = float(line.split(':')[1].strip())
        else:
            print(f"Warning: lr_leaveout.log not found after testing without '{f}'")

        # Store results
        df_r2.loc[f, i] = r2_val
        df_mae.loc[f, i] = mae_val
        df_mse.loc[f, i] = mse_val

    # Choose the next feature to remove based on current round (e.g., highest MAE)
    next_feature = df_mae[i].astype(float).idxmax()
    print(f"Iteration {i}: removing '{next_feature}'")
    remaining_features.remove(next_feature)

# Save results
df_r2.to_csv('leaveout_r2.tsv', sep='\t')
df_mae.to_csv('leaveout_mae.tsv', sep='\t')
df_mse.to_csv('leaveout_mse.tsv', sep='\t')

# import pandas as pd
# import matplotlib.pyplot as plt

# # Load results
# df_mae = pd.read_csv('leaveout_mae.tsv', sep='\t', index_col=0)
# df_mse = pd.read_csv('leaveout_mse.tsv', sep='\t', index_col=0)

# Get the feature removal order from df_mae columns (i.e., by iteration)
removal_order = []
for col in df_mae.columns:
    col_series = df_mae[col].dropna()
    if not col_series.empty:
        # The removed feature this iteration is the row with the highest MAE
        removed = col_series.astype(float).idxmax()
        removal_order.append(removed)

# Collect the corresponding MAE and MSE values from the selected removal features
mae_values = [df_mae.loc[feat].dropna().iloc[-1] for feat in removal_order]
mse_values = [df_mse.loc[feat].dropna().iloc[-1] for feat in removal_order]

# Plot MAE
plt.figure(figsize=(12, 5))
plt.plot(range(len(mae_values)), mae_values, marker='o', label='MAE')
plt.xticks(range(len(removal_order)), removal_order, rotation=90)
plt.ylabel('MAE')
plt.xlabel('Removed Features (by iteration)')
plt.title('MAE vs Feature Removal Order')
plt.tight_layout()
plt.savefig('leaveout_mae.png')
plt.close()

# Plot MSE
plt.figure(figsize=(12, 5))
plt.plot(range(len(mse_values)), mse_values, marker='o', label='MSE', color='orange')
plt.xticks(range(len(removal_order)), removal_order, rotation=90)
plt.ylabel('MSE')
plt.xlabel('Removed Features (by iteration)')
plt.title('MSE vs Feature Removal Order')
plt.tight_layout()
plt.savefig('leaveout_mse.png')
plt.close()
