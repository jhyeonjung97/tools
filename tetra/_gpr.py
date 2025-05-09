import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ylabels = {
    'coord': 'Coordination',
    'row': 'Row',
    'numb': 'Number',
    'metal': 'Metal',
    'CN': 'Coordination Number',
    'ON': 'Oxidation Number',
    'energy': 'Energy (eV)',
    'form': 'Formation Energy (eV)',
    'coh': 'Cohesive Energy (eV)',
    'volume': 'Volume (Å³)',
    'cell': 'Cell',
    'chg': 'Bader Charge (e⁻)',
    'mag': 'Magnetic Moments (μB)',
    'l_bond': 'Bond Length (Å)',
    'n_bond': 'Number of Bonds per Metal',
    'match': 'Bulk Structure Maintain',
    '-ICOHP': '-ICOHP per Metal (eV)',
    'ICOBI': 'ICOBI per Metal',
    '-ICOOP': '-ICOOP per Metal (eV)',
    '-ICOHPn': '-ICOHP per Bond (eV)',
    'ICOBIc': 'ICOBI per Bond',
    '-ICOOPc': '-ICOOP per Bond (eV)',
    'madelung': 'Madelung Energy (Loewdin eV)',
    'grosspop': 'Gross Population (Loewdin e⁻)',
}

def main():
    parser = argparse.ArgumentParser(description='Gaussian Process regression using bulk_data.csv and mendeleev_data.csv')
    parser.add_argument('--Y', default='form', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=[
        'chg', 'mag', 'volume', 'l_bond', 'n_bond',
        'grosspop', 'madelung', 'ICOHP', 'ICOHPc', 'ICOBI', 'ICOBIc', 'ICOOP', 'ICOOPc', 
        'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
        'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    args = parser.parse_args()

    # Convert feature names if they start with ICOHP or ICOOP (prepend '-')
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]

    root = '/Users/hailey/Desktop/7_V_bulk/figures'
    df_bulk = pd.read_csv(os.path.join(root, 'bulk_data.csv'), index_col=0)
    df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'), index_col=0)

    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    df = df.rename(columns={'row_bulk': 'row', 'numb_mend': 'numb'})
    df = df.drop(columns=['row_mend', 'numb_bulk'])
    df = df[df['row'] != 'fm']

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    all_columns = args.X + [args.Y]
    df = df.dropna(subset=all_columns)

    X = df[args.X].astype(float).values
    Y = df[args.Y].astype(float).values

    kernel = RBF(length_scale=1.0) + WhiteKernel(noise_level=1.0)
    model = GaussianProcessRegressor(kernel=kernel, normalize_y=True, random_state=42)

    model.fit(X, Y)
    Y_pred, Y_std = model.predict(X, return_std=True)

    mae = mean_absolute_error(Y, Y_pred)
    mse = mean_squared_error(Y, Y_pred)
    r2 = r2_score(Y, Y_pred)

    output_suffix = args.output

    # Save regression summary
    with open(f'gp_{output_suffix}.log', 'w') as f:
        f.write("Gaussian Process Regressor\n")
        f.write(f"Kernel: {model.kernel_}\n")
        f.write(f"R2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['Y_std'] = Y_std
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(f'gp_{output_suffix}.tsv', sep='\t', index=False)

    # Plot parity
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}
    coord_map = {'WZ': '>', 'ZB': '<', 'TN': 'o', 'PD': 'o', 'NB': 's', 'RS': 'd', 'LT': 'h'}

    plt.figure(figsize=(10, 8))
    for r in df['row'].unique():
        for c in df['coord'].unique():
            subset = df[(df['row'] == r) & (df['coord'] == c)]
            if subset.empty:
                continue
            X_subset = subset[args.X].astype(float).values
            Y_subset = subset[args.Y].values
            Y_pred_subset = model.predict(X_subset)
            plt.scatter(
                Y_subset,
                Y_pred_subset,
                label=f'{r}_{c}',
                alpha=0.3,
                color=row_map.get(r, 'gray'),
                marker=coord_map.get(c, 'x')
            )
            for _, row_data in subset.iterrows():
                x_single = np.array([row_data[args.X].values.astype(float)])
                y_pred_single = model.predict(x_single)[0]
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single), fontsize=8)

    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels.get(args.Y, args.Y)}')
    plt.ylabel(f'Predicted {ylabels.get(args.Y, args.Y)}')
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(f'gp_{output_suffix}.png')
    plt.close()

    print(f"Saved: gp_{output_suffix}.log, gp_{output_suffix}.tsv, gp_{output_suffix}.png")

if __name__ == '__main__':
    main()
