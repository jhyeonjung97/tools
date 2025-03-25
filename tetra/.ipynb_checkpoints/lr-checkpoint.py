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
        'CN', 'ON', 'energy', 'volume', 'cell', 'chg', 'mag', 'l_bond', 'n_bond',
        '-ICOHPm', 'ICOBIm', '-ICOOPm', '-ICOHPn', 'ICOBIn', '-ICOOPn',
        'madelung', 'grosspop', 'Natom', 'Vatom', 'Tboil', 'Tmelt', 'mass', 'density',
        'dipole', 'pauling', 'Rcoval', 'Rmetal', 'Rvdw', 'Hevap', 'Hfus', 'Hform',
        'ion1', 'ion2', 'ion3', 'ion12'
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename prefix')
    args = parser.parse_args()

    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
    df_bulk = pd.read_csv(os.path.join(root, 'bulk_data.csv'), index_col=0)
    df_bulk['metal'] = df_bulk['metal'].astype(str)

    # Keep row and coord as string for filtering
    df_bulk['row'] = df_bulk['row'].astype(str)
    df_bulk['coord'] = df_bulk['coord'].astype(str)

    df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'), index_col=0)

    # Merge
    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    
    print(df.head(10))

    # Optional filtering
    if args.row:
        df = df[df['row'] == args.row]
    if args.coord:
        df = df[df['coord'] == args.coord]

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

    output_prefix = args.output

    # Save regression summary
    with open(f'{output_prefix}.log', 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(f'{output_prefix}.tsv', sep='\t', index=False)

    # Plot parity
    plt.figure(figsize=(6, 6))
    plt.scatter(Y, Y_pred, alpha=0.7)
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], 'r--')
    plt.xlabel('DFT-calculated')
    plt.ylabel('Predicted')
    plt.title(f'{args.Y} prediction')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png')
    plt.close()

    print(f"Saved: {output_prefix}.log, .tsv, .png")


if __name__ == '__main__':
    main()

    parser.add_argument('--row', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename prefix')
    args = parser.parse_args()

    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
    df_bulk = pd.read_csv(os.path.join(root, 'bulk_data.csv'))
    df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'))

    # Merge on metal
    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    
    print(df)

    # Optional filtering
    if args.row:
        df = df[df['row'] == args.row]
    if args.coord:
        df = df[df['coord'] == args.coord]

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

    output_prefix = args.output

    # Save regression summary
    with open(f'{output_prefix}.log', 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(f'{output_prefix}.tsv', sep='\t', index=False)

    # Plot parity
    plt.figure(figsize=(6, 6))
    plt.scatter(Y, Y_pred, alpha=0.7)
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], 'r--')
    plt.xlabel('DFT-calculated')
    plt.ylabel('Predicted')
    plt.title(f'{args.Y} prediction')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png')
    plt.close()

    print(f"Saved: {output_prefix}.log, .tsv, .png")


if __name__ == '__main__':
    main()
