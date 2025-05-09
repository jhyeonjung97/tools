import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import socket

hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk/figures'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
elif user_name == 'hailey':
    root = '/Users/hailey/Desktop/7_V_bulk/figures'
else:
    raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

# Define the full feature list
features = [
    'numb', 'chg', 'mag', 'volume', 'l_bond', 'n_bond', 'grosspop', 'madelung', 
    'ICOHP', 'ICOHPc', 'ICOBI', 'ICOBIc', 'ICOOP', 'ICOOPc', 
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

for i in range(len(features)-1):
    for f in remaining_features:
        test_features = selected_features + [f]
        input_str = ' '.join(test_features)

        # Remove previous log if exists
        log_path = os.path.join(root, 'lr_addon.log')
        if os.path.exists(log_path):
            os.remove(log_path)

        # Run lr.py with selected features
        cmd = f'python ~/bin/tools/tetra/lr.py --Y form --X {input_str} --output addon'
        subprocess.run(cmd, shell=True, check=True)

        # Parse lr_addon.log
        r2_val = mae_val = mse_val = None
        if os.path.exists(log_path):
            with open(log_path, 'r') as log_file:
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
df_r2.to_csv(os.path.join(root, 'addon_r2.tsv'), sep='\t')
df_mae.to_csv(os.path.join(root, 'addon_mae.tsv'), sep='\t')
df_mse.to_csv(os.path.join(root, 'addon_mse.tsv'), sep='\t')

# import pandas as pd
# import matplotlib.pyplot as plt

# # Load results
# df_mae = pd.read_csv('addon_mae.tsv', sep='\t', index_col=0)
# df_mse = pd.read_csv('addon_mse.tsv', sep='\t', index_col=0)

# Find order of selected features from df_mae (column-wise non-NaN values)
selected_features_order = []
for col in df_mae.columns:
    col_series = df_mae[col].dropna()
    if not col_series.empty:
        # The selected feature this iteration is the row with smallest MAE
        selected = col_series.astype(float).idxmin()
        selected_features_order.append(selected)

# Gather MAE and MSE values in order of feature addition
mae_values = [df_mae.loc[feat].dropna().iloc[-1] for feat in selected_features_order]
mse_values = [df_mse.loc[feat].dropna().iloc[-1] for feat in selected_features_order]

# Plot MAE
plt.figure(figsize=(12, 5))
plt.plot(range(len(mae_values)), mae_values, marker='o', label='MAE')
plt.xticks(range(len(selected_features_order)), selected_features_order, rotation=90)
plt.ylabel('MAE')
plt.xlabel('Added Features (by iteration)')
plt.title('MAE vs Feature Addition Order')
plt.tight_layout()
plt.savefig(os.path.join(root, 'addon_mae.png'))
plt.close()

# Plot MSE
plt.figure(figsize=(12, 5))
plt.plot(range(len(mse_values)), mse_values, marker='o', label='MSE', color='orange')
plt.xticks(range(len(selected_features_order)), selected_features_order, rotation=90)
plt.ylabel('MSE')
plt.xlabel('Added Features (by iteration)')
plt.title('MSE vs Feature Addition Order')
plt.tight_layout()
plt.savefig(os.path.join(root, 'addon_mse.png'))
plt.close()
