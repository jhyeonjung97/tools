import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
import socket

hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk/figures'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/7_V_bulk/figures'
else:
    raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

ylabels = {
    'form': 'Formation Energy (eV)',
    'coh': 'Cohesive Energy (eV)',
}

def select_features_by_correlation(correlation_matrix, target_col, threshold=0.7):
    """
    Select one feature from a group of highly correlated features
    
    Args:
        correlation_matrix: Correlation matrix of features
        target_col: Name of the target variable
        threshold: Threshold for high correlation (default: 0.7)
    
    Returns:
        selected_features: List of selected features
    """
    # Sort features by absolute correlation with target variable
    target_corr = correlation_matrix[target_col].abs().sort_values(ascending=False)
    
    # Track selected and dropped features
    selected_features = []
    dropped_features = []
    
    # Iterate through all features
    for feature in target_corr.index:
        if feature == target_col:
            continue
            
        # Check correlation with already selected features
        high_corr = False
        for selected in selected_features:
            if abs(correlation_matrix.loc[feature, selected]) > threshold:
                high_corr = True
                # Select feature with higher correlation to target
                if target_corr[feature] > target_corr[selected]:
                    selected_features.remove(selected)
                    dropped_features.append(selected)
                    selected_features.append(feature)
                else:
                    dropped_features.append(feature)
                break
                
        if not high_corr:
            selected_features.append(feature)
    
    print("\nSelected features:")
    for feat in selected_features:
        print(f"- {feat}")
    
    print("\nDropped features:")
    for feat in dropped_features:
        print(f"- {feat}")
    
    return selected_features

def main():
    parser = argparse.ArgumentParser(description='Linear regression using bulk_data.csv and mendeleev_data.csv')
    parser.add_argument('--Y', default='form', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=[
        'OS', 'CN', 'numb', 'chg', 'mag', 'volume', 'l_bond', 'madelung',
        'ICOHPm', 'ICOHPmn', 'ICOHPn', 'ICOBIm', 'ICOBImn', 'ICOBIn', 'ICOOPm', 'ICOOPmn', 'ICOOPn', 
        # 'ICOHPm', 'ICOHPn', 'ICOBIm', 'ICOBIn', 'ICOOPm', 'ICOOPn', 
        'ion-1', 'ion', 'ion+1', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1', 
        'pauling', 'Natom', 'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    args = parser.parse_args()
    
    # Convert feature names if they start with ICOHP or ICOOP (prepend '-')
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    
    # Read data files
    df_tetra = pd.read_csv(os.path.join(root, 'bulk_data.csv'), index_col=0, dtype={'coord': str})
    df_octa = pd.read_csv(os.path.join(root, 'comer_bulk_data.csv'), index_col=0, dtype={'coord': str})
    df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'), index_col=0)

    df_bulk = pd.concat([df_tetra, df_octa], axis=0)
    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    df = df.rename(columns={'row_bulk': 'row', 'numb_mend': 'numb'})
    df = df.drop(columns=['row_mend', 'numb_bulk'])
    df = df[df['row'] != 'fm'] ##

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    # Calculate ion values for each row
    df['ion-1'] = df.apply(lambda row: row[f'ion{int(row["OS"])-1}'], axis=1)
    df['ion'] = df.apply(lambda row: row[f'ion{int(row["OS"])}'], axis=1)
    df['ion+1'] = df.apply(lambda row: row[f'ion{int(row["OS"])+1}'], axis=1)

    df['ion-1n'] = df['ion-1']/df['OS']
    df['ionn'] = df['ion']/df['OS']
    df['ion+1n'] = df['ion+1']/df['OS']

    df['-ICOHPmn'] = df['-ICOHPm']/df['OS']
    df['ICOBImn'] = df['ICOBIm']/df['OS']
    df['-ICOOPmn'] = df['-ICOOPm']/df['OS']

    # Initialize columns with float type
    df['ionN-1'] = 0.0
    df['ionN'] = 0.0
    df['ionN+1'] = 0.0
    
    for idx, row in df.iterrows():
        for i in range(1, int(row['OS'])):
            df.at[idx, 'ionN-1'] = float(df.at[idx, 'ionN-1']) + float(row[f'ion{i}'])
        df.at[idx, 'ionN'] = float(df.at[idx, 'ionN-1']) + float(row[f'ion{int(row["OS"])}'])
        df.at[idx, 'ionN+1'] = float(df.at[idx, 'ionN']) + float(row[f'ion{int(row["OS"])+1}'])

    # Drop rows with NaN in any relevant column
    df = df.dropna(subset=args.X + [args.Y])

    # Exclude ion1~7 columns from saving
    columns_to_save = [col for col in ['metal', 'row', 'coord'] + args.X + [args.Y] if col not in ['ion1', 'ion2', 'ion3', 'ion4', 'ion5', 'ion6', 'ion7']]
    df_to_save = df[columns_to_save]
    df_to_save.to_csv(f'{root}/bulk_data_total.csv', sep=',')
    df_to_save.to_csv(f'{root}/bulk_data_total.tsv', sep='\t', float_format='%.2f')

    X = df[args.X].astype(float)
    Y = df[args.Y].astype(float)

    model = LinearRegression()
    model.fit(X, Y)

    Y_pred = model.predict(X)
    mae = mean_absolute_error(Y, Y_pred)
    mse = mean_squared_error(Y, Y_pred)
    r2 = model.score(X, Y)

    output_suffix = args.output

    # Save regression summary
    with open(os.path.join(root, f'lr_{output_suffix}.log'), 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(os.path.join(root, f'lr_{output_suffix}.tsv'), sep='\t', index=False)

    # Plot parity with color by 'row' and marker by 'coord'
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}
    coord_map = {'WZ': '+', 'ZB': 'x', 'TN': 'o', 'PD': 'o', 'NB': 's', 'RS': 'D', 'LT': 'h', '+3': 'v', '+4': '^', '+5': '<', '+6': '>'}

    plt.figure(figsize=(10, 8))
    for r in df['row'].unique():
        for c in df['coord'].unique():
            subset = df[(df['row'] == r) & (df['coord'] == c)]
            if subset.empty:
                continue
            plt.scatter(
                subset[args.Y],
                model.predict(subset[args.X].astype(float)),
                label=f'{r}_{c}',
                alpha=0.3,
                color=row_map.get(r, 'gray'),
                marker=coord_map.get(c, 'x')
            )
            for _, row_data in subset.iterrows():
                row_features = pd.DataFrame([row_data[args.X].values], columns=args.X)
                y_pred_single = model.predict(row_features)[0]
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single), fontsize=8)

    
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(root, f'lr_{output_suffix}.png'))
    plt.close()

    # Save covariance and correlation matrices
    df_metrics = pd.concat([Y, X], axis=1)  # Place Y first
    cov = df_metrics.cov()
    cor = df_metrics.corr()

    # Save results using all features
    print("\nTesting with all features")
    model_all = LinearRegression()
    model_all.fit(X, Y)
    
    Y_pred_all = model_all.predict(X)
    mae_all = mean_absolute_error(Y, Y_pred_all)
    mse_all = mean_squared_error(Y, Y_pred_all)
    r2_all = model_all.score(X, Y)
    
    # Save results for all features
    with open(os.path.join(root, f'lr_{args.output}_all.log'), 'w') as f:
        f.write(f"Using all features\n")
        f.write(f"Intercept: {model_all.intercept_:.4f}\n")
        for name, coef in zip(args.X, model_all.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2_all:.4f}\nMAE: {mae_all:.4f}\nMSE: {mse_all:.4f}\n")

    # Save data with predictions
    df_result_all = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result_all['Y_true'] = Y
    df_result_all['Y_pred'] = Y_pred_all
    df_result_all['residual'] = Y - Y_pred_all
    df_result_all.to_csv(os.path.join(root, f'lr_{args.output}_all.tsv'), sep='\t', index=False)

    plt.figure(figsize=(10, 8))
    for r in df['row'].unique():
        for c in df['coord'].unique():
            subset = df[(df['row'] == r) & (df['coord'] == c)]
            if subset.empty:
                continue
            plt.scatter(
                subset[args.Y],
                model_all.predict(subset[args.X].astype(float)),
                label=f'{r}_{c}',
                alpha=0.3,
                color=row_map.get(r, 'gray'),
                marker=coord_map.get(c, 'x')
            )
            for _, row_data in subset.iterrows():
                row_features = pd.DataFrame([row_data[args.X].values], columns=args.X)
                y_pred_single = model_all.predict(row_features)[0]
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single), fontsize=8)

    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(root, f'lr_{args.output}_all.png'))
    plt.close()

    # Save covariance and correlation matrices
    cov.to_csv(os.path.join(root, f'covariance_{args.output}_all.tsv'), sep='\t')
    cor.to_csv(os.path.join(root, f'correlation_{args.output}_all.tsv'), sep='\t')

    plt.figure(figsize=(10, 8))
    sns.heatmap(cov, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 6},
                cbar_kws={"shrink": 0.5, "aspect": 20})  # Adjust colorbar width
    plt.tight_layout(pad=1.0)  # Reduce padding
    plt.savefig(os.path.join(root, f'covariance_{args.output}_all.png'))
    plt.close()

    plt.figure(figsize=(10, 8))
    sns.heatmap(cor, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 6},
                cbar_kws={"shrink": 0.5, "aspect": 20})  # Adjust colorbar width
    plt.tight_layout(pad=1.0)  # Reduce padding
    plt.savefig(os.path.join(root, f'correlation_{args.output}_all.png'))
    plt.close()

    print(f"Saved results using all features")

    # Test different correlation thresholds
    thresholds = [0.7, 0.8, 0.9]
    threshold_suffixes = ['7', '8', '9']

    for threshold, suffix in zip(thresholds, threshold_suffixes):
        print(f"\nTesting correlation threshold: {threshold}")
        
        # Select features based on correlation threshold
        selected_features = select_features_by_correlation(cor, args.Y, threshold)
        
        # Train new model with selected features
        X_selected = df[selected_features].astype(float)
        model_selected = LinearRegression()
        model_selected.fit(X_selected, Y)
        
        Y_pred_selected = model_selected.predict(X_selected)
        mae_selected = mean_absolute_error(Y, Y_pred_selected)
        mse_selected = mean_squared_error(Y, Y_pred_selected)
        r2_selected = model_selected.score(X_selected, Y)
        
        # Save results for selected features
        with open(os.path.join(root, f'lr_{args.output}_{suffix}.log'), 'w') as f:
            f.write(f"Correlation threshold: {threshold}\n")
            f.write(f"Selected Features: {', '.join(selected_features)}\n")
            f.write(f"Intercept: {model_selected.intercept_:.4f}\n")
            for name, coef in zip(selected_features, model_selected.coef_):
                f.write(f"{name}: {coef:.4f}\n")
            f.write(f"\nR2: {r2_selected:.4f}\nMAE: {mae_selected:.4f}\nMSE: {mse_selected:.4f}\n")

        # Save data with predictions
        df_result_selected = df[['metal', 'row', 'coord'] + selected_features].copy()
        df_result_selected['Y_true'] = Y
        df_result_selected['Y_pred'] = Y_pred_selected
        df_result_selected['residual'] = Y - Y_pred_selected
        df_result_selected.to_csv(os.path.join(root, f'lr_{args.output}_{suffix}.tsv'), sep='\t', index=False)

        plt.figure(figsize=(10, 8))
        for r in df['row'].unique():
            for c in df['coord'].unique():
                subset = df[(df['row'] == r) & (df['coord'] == c)]
                if subset.empty:
                    continue
                plt.scatter(
                    subset[args.Y],
                    model_selected.predict(subset[selected_features].astype(float)),
                    label=f'{r}_{c}',
                    alpha=0.3,
                    color=row_map.get(r, 'gray'),
                    marker=coord_map.get(c, 'x')
                )
                for _, row_data in subset.iterrows():
                    row_features = pd.DataFrame([row_data[selected_features].values], columns=selected_features)
                    y_pred_single = model_selected.predict(row_features)[0]
                    plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single), fontsize=8)

        plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
        plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
        plt.ylabel(f'Predicted {ylabels[args.Y]}')
        plt.legend(loc='best', fontsize=8)
        plt.tight_layout()
        plt.savefig(os.path.join(root, f'lr_{args.output}_{suffix}.png'))
        plt.close()

        # Save covariance and correlation matrices for selected features
        df_metrics_selected = pd.concat([Y, X_selected], axis=1)  # Place Y first
        cov_selected = df_metrics_selected.cov()
        cor_selected = df_metrics_selected.corr()

        cov_selected.to_csv(os.path.join(root, f'covariance_{args.output}_{suffix}.tsv'), sep='\t')
        cor_selected.to_csv(os.path.join(root, f'correlation_{args.output}_{suffix}.tsv'), sep='\t')

        plt.figure(figsize=(10, 8))
        sns.heatmap(cov_selected, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 6},
                    cbar_kws={"shrink": 0.5, "aspect": 20})  # Adjust colorbar width
        plt.tight_layout(pad=1.0)  # Reduce padding
        plt.savefig(os.path.join(root, f'covariance_{args.output}_{suffix}.png'))
        plt.close()

        plt.figure(figsize=(10, 8))
        sns.heatmap(cor_selected, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 6},
                    cbar_kws={"shrink": 0.5, "aspect": 20})  # Adjust colorbar width
        plt.tight_layout(pad=1.0)  # Reduce padding
        plt.savefig(os.path.join(root, f'correlation_{args.output}_{suffix}.png'))
        plt.close()

        print(f"Saved results for threshold {threshold} with suffix {suffix}")

    print(f"Saved: lr_{output_suffix}.log, lr_{output_suffix}.tsv")
    print(f"Saved: covariance_{output_suffix}.tsv, covariance_{output_suffix}.png")
    print(f"Saved: correlation_{output_suffix}.tsv, correlation_{output_suffix}.png")


if __name__ == '__main__':
    main()