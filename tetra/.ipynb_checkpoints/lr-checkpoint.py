import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
import os


def main():
    parser = argparse.ArgumentParser(description='Linear regression using bulk_data.csv and mendeleev_data.csv')
    parser.add_argument('--Y', default='form', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=[
        'chg', 'mag', 'volume', 'l_bond', 'n_bond',
        'grosspop', 'madelung', '-ICOHPm', '-ICOHPn', 'ICOBIm', 'ICOBIn', '-ICOOPm', '-ICOOPn', 
        'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
        'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    args = parser.parse_args()

    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
    df_bulk = pd.read_csv(os.path.join(root, 'bulk_data.csv'), index_col=0)
    df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'), index_col=0)

    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    df = df.rename(columns={'row_bulk': 'row', 'numb_mend': 'numb'})
    df = df.drop(columns=['row_mend', 'numb_bulk'])
    df = df[df['coord'] != 'fm']

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    # Drop rows with NaN in any relevant column
    all_columns = args.X + [args.Y]
    df = df.dropna(subset=all_columns)

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
    with open(f'lr_{output_suffix}.log', 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(f'lr_{output_suffix}.tsv', sep='\t', index=False)

    # Plot parity with color by 'row' and marker by 'coord'
    colors = ['red', 'green', 'blue']
    markers = ['>', '<', 'o', 's', 'p', 'd', 'h', '^', 'v']
    row_map = {r: colors[i] for i, r in enumerate(sorted(df['row'].unique()))}
    coord_map = {c: markers[i] for i, c in enumerate(sorted(df['coord'].unique()))}

    plt.figure(figsize=(6, 6))
    for r in df['row'].unique():
        for c in df['coord'].unique():
            subset = df[(df['row'] == r) & (df['coord'] == c)]
            plt.scatter(
                subset[args.Y],
                model.predict(subset[args.X].astype(float)),
                label=f'{r}_{c}',
                alpha=0.7,
                color=row_map[r],
                marker=coord_map[c]
            )

    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], 'r--')
    plt.xlabel('DFT-calculated')
    plt.ylabel('Predicted')
    plt.title(f'{args.Y} prediction')
    plt.legend(loc='best', fontsize=6)
    plt.tight_layout()
    plt.savefig(f'lr_{output_suffix}.png')
    plt.close()

    # Save covariance and correlation matrices as TSV and PNG
    df_metrics = pd.concat([Y, X], axis=1)  # Place Y first
    cov = df_metrics.cov()
    cor = df_metrics.corr()

    cov.to_csv(f'covariance_{output_suffix}.tsv', sep='\t')
    cor.to_csv(f'correlation_{output_suffix}.tsv', sep='\t')

    plt.figure(figsize=(10, 8))
    sns.heatmap(cov, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 5})
    plt.tight_layout()
    plt.savefig(f'covariance_{output_suffix}.png')
    plt.close()

    plt.figure(figsize=(10, 8))
    sns.heatmap(cor, annot=True, fmt='.2f', cmap='coolwarm', annot_kws={"size": 5})
    plt.tight_layout()
    plt.savefig(f'correlation_{output_suffix}.png')
    plt.close()

    print(f"Saved: {output_suffix}.log, {output_suffix}.tsv, {output_suffix}.png")
    print(f"Saved: {output_suffix}_covariance.tsv, {output_suffix}_covariance.png")
    print(f"Saved: {output_suffix}_correlation.tsv, {output_suffix}_correlation.png")


if __name__ == '__main__':
    main()