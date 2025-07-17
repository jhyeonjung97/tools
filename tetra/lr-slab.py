import os
import socket
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.linear_model import LinearRegression
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import mean_absolute_error, mean_squared_error
from mendeleev import element
import math

def main():
    hostname = socket.gethostname()
    user_name = os.getlogin()
    if hostname == 'PC102616':
        root = '/Users/jiuy97/Desktop'
    elif user_name == 'jiuy97':
        root = '/pscratch/sd/j/jiuy97'
    elif user_name == 'hailey' or user_name == 'root':
        root = '/Users/hailey/Desktop'
    else:
        raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

    bulk_dir = os.path.join(root, '7_V_bulk/figures')
    slab_dir = os.path.join(root, '8_V_slab/figures')

    coords_data = [
        {'coord': 'WZ', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.667, 'coord_dir': '1_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
        {'coord': 'ZB', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.833, 'coord_dir': '2_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
        {'coord': 'TN', 'CN': 4, 'OS': 2, 'MN': 4, 'surfox': 1.933, 'coord_dir': '3_SquarePlanar_TN', 'zorder': 3, 'marker': 'o', 'color': 'dodgerblue',},
        {'coord': 'PD', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.667, 'coord_dir': '4_SquarePlanar_PD', 'zorder': 2, 'marker': 'o', 'color': 'deepskyblue',},
        {'coord': 'NB', 'CN': 4, 'OS': 2, 'MN': 6, 'surfox': 2.667, 'coord_dir': '5_SquarePlanar_NB', 'zorder': 1, 'marker': 's', 'color': 'limegreen',},
        {'coord': 'RS', 'CN': 6, 'OS': 2, 'MN': 2, 'surfox': 1.933, 'coord_dir': '6_Octahedral_RS',   'zorder': 6, 'marker': 'd', 'color': 'orchid',},
        {'coord': 'LT', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '7_Pyramidal_LT',    'zorder': 0, 'marker': 'h', 'color': 'silver',},
        {'coord': 'WW', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '8_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
        {'coord': 'ZZ', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '9_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
        {'coord': '+3', 'CN': 6, 'OS': 3, 'MN': 4, 'surfox': 2.833, 'coord_dir': '1_Octahedral_+3_012', 'zorder': 1, 'marker': 'o', 'color': 'darkorange'},
        {'coord': '+4', 'CN': 6, 'OS': 4, 'MN': 2, 'surfox': 4.000, 'coord_dir': '2_Octahedral_+4_100', 'zorder': 2, 'marker': 's', 'color': 'gold'},
        {'coord': '+4', 'CN': 6, 'OS': 4, 'MN': 2, 'surfox': 3.333, 'coord_dir': '3_Octahedral_+4_110', 'zorder': 2, 'marker': 's', 'color': 'silver'},
        {'coord': '+5', 'CN': 6, 'OS': 5, 'MN': 4, 'surfox': 5.000, 'coord_dir': '4_Octahedral_+5_100', 'zorder': 3, 'marker': '^', 'color': 'dodgerblue'},
        {'coord': '+6', 'CN': 6, 'OS': 6, 'MN': 1, 'surfox': 6.000, 'coord_dir': '5_Octahedral_+6_001', 'zorder': 4, 'marker': 'v', 'color': 'deepskyblue'},
    ]
    coords = pd.DataFrame(coords_data).set_index('coord_dir')
    coords.index.name = None
    coord_order = [item['coord'] for item in coords_data]
    coord_order_dict = {coord: idx for idx, coord in enumerate(coord_order)}

    metals = {
        '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
        '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
        '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    }
    indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])]

    X_default = [
        'bulkox', 'bulk_e', 'bulk_ie', 'bulk_ieo', 'bulk_ien', 'bulk_cfse',
        'surfox', 'surf_e', 'surf_ie', 'surf_ieo', 'surf_ien', 'surf_cfse',
        'ICOHP', 'ICOHPo', 'ICOHPn', 'Hform', 'Hsub', 'Hevap', 'Hfus',
        'CN', 'group', 'period', 'chg', 'mag', 'volume', 'l_bond', 'cfse',
    ]
    
    ylabels = {
        'o_energy': 'O Adsorption Energy (eV)',
        'oh_energy': 'OH Adsorption Energy (eV)',
    }

    parser = argparse.ArgumentParser()
    parser.add_argument('--Y', default='o_energy', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=X_default, help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    args = parser.parse_args()

    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]

    df_csv = os.path.join(slab_dir, 'lr_slab_data.csv')
    df_tsv = os.path.join(slab_dir, 'lr_slab_data.tsv')
    if os.path.exists(df_csv):
        df = pd.read_csv(df_csv, index_col=0)
    else:
        df_mend = pd.read_csv(os.path.join(bulk_dir, 'mendeleev_data.csv'), index_col=0)
        df_bulk_tetra = pd.read_csv(os.path.join(bulk_dir, 'bulk_data.csv'), index_col=0, dtype={'coord': str})
        df_bulk_octa = pd.read_csv(os.path.join(bulk_dir, 'comer_bulk_data.csv'), index_col=0, dtype={'coord': str})
        df_slab_tetra = pd.read_csv(os.path.join(slab_dir, 'slab_data.csv'), index_col=0, dtype={'coord': str})
        df_slab_octa = pd.read_csv(os.path.join(slab_dir, 'slab_data_comer.csv'), index_col=0, dtype={'coord': str})

        df_bulk = pd.concat([df_bulk_tetra, df_bulk_octa], axis=0)
        df_slab = pd.concat([df_slab_tetra, df_slab_octa], axis=0)
        df_bulk = df_bulk.sort_values(['coord', 'row', 'numb'], key=lambda x: x.map(coord_order_dict) if x.name == 'coord' else x)
        df_slab = df_slab.sort_values(['coord', 'row', 'numb'], key=lambda x: x.map(coord_order_dict) if x.name == 'coord' else x)

        df_sum = pd.merge(df_slab, df_bulk, left_index=True, right_index=True, suffixes=('_slab', '_bulk'))
        df_sum = df_sum.rename(columns={'coord_slab': 'coord', 'row_slab': 'row', 'numb_slab': 'numb', 'metal_slab': 'metal'})
        df_sum = df_sum.drop(columns=['coord_bulk', 'row_bulk', 'numb_bulk', 'metal_bulk'])

        df = pd.merge(df_sum, df_mend, left_on='metal', right_index=True, suffixes=('_slab', '_mend'))
        df = df.rename(columns={'row_slab': 'row', 'numb_slab': 'numb', 'group_slab': 'group'})
        df = df.drop(columns=['row_mend', 'numb_mend', 'group_mend'])

        row2period = {
            '3d': 4,
            '4d': 5,
            '5d': 6,
        }
        df['period'] = df['row'].map(row2period)
        df['Hsub'] = df['Hevap'] + df['Hfus']

        df['chgc'] = df['chg'] / df['CN']
        df['chgo'] = df['chg'] / df['OS']
        df['chgn'] = df['chg'] / df['OS'] / df['CN']

        df['-ICOHPo'] = df['-ICOHP'] / df['OS']
        df['ICOBIo'] = df['ICOBI'] / df['OS']
        df['-ICOOPo'] = df['-ICOOP'] / df['OS']

        df['-ICOHPn'] = df['-ICOHPo'] / df['CN']
        df['ICOBIn'] = df['ICOBIo'] / df['CN']
        df['-ICOOPn'] = df['-ICOOPo'] / df['CN']

        df['bulkox'] = df['OS']
        df['surfox_d'] = df['surfox'] - np.floor(df['surfox'])

        df['bulk_e'] = df['group'] - df['bulkox']
        df['surf_e'] = df['group'] - df['surfox']
        df['surf_ef'] = df.apply(lambda row: math.floor(row['surf_e']+1), axis=1)
        df['surf_ec'] = df.apply(lambda row: math.ceil(row['surf_e']+1), axis=1)
        df = df[(df['bulk_e'] >= 0) & (df['surf_e'] >= 0)]

        df['ions'] = df['metal'].apply(lambda x: element(x).ionenergies)
        df['bulk_ie'] = df.apply(lambda row: row['ions'].get(row['bulk_e']+1, 0) / row['bulk_e'] if row['bulk_e'] != 0 else 0, axis=1)
        df['surf_ie'] = df.apply(lambda row: ((row['ions'].get(row['surf_ef'], 0) * (1-row['surfox_d']) + row['ions'].get(row['surf_ec'], 0) * row['surfox_d']) / row['surfox']), axis=1)
        df['bulk_ieo'] = df['bulk_ie'] / df['OS']
        df['surf_ieo'] = df['surf_ie'] / df['OS']
        df['bulk_ien'] = df['bulk_ie'] / df['OS'] / df['CN']
        df['surf_ien'] = df['surf_ie'] / df['OS'] / df['CN']

        df['bulk_de'] = df['bulk_e'].apply(lambda x: max(0, min(10, round(x))))
        df['surf_de'] = df['surf_e'].apply(lambda x: max(0, min(10, round(x))))
        df['bulk_cfse'] = df.apply(lambda row: calculate_cfse(row['coord'], row['bulk_de'], row['row']), axis=1)
        df['surf_cfse'] = df.apply(lambda row: calculate_cfse(row['coord'], row['surf_de'], row['row']), axis=1)
        df['cfse'] = df['bulk_cfse'] - df['surf_cfse']

        df = df.drop(columns=[col for col in df.columns if col not in args.X + ['metal', 'row', 'coord', 'o_energy', 'oh_energy']])

        df.to_csv(df_csv)
        df.to_csv(df_tsv, sep='\t', float_format='%.2f')

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    df = df.dropna(subset=args.X + [args.Y])
    print(df)

    X = df[args.X].astype(float)
    Y = df[args.Y].astype(float)

    model = LinearRegression()
    model.fit(X, Y)

    Y_pred = model.predict(X)
    mae = mean_absolute_error(Y, Y_pred)
    mse = mean_squared_error(Y, Y_pred)
    r2 = model.score(X, Y)

    if args.Y == 'o_energy':
        output_suffix = 'o_' + args.output
    elif args.Y == 'oh_energy':
        output_suffix = 'oh_' + args.output
    else:
        output_suffix = args.output

    # Save regression summary
    with open(os.path.join(slab_dir, f'lr_slab_{output_suffix}.log'), 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")
    print(f"saved log in {os.path.join(slab_dir, f'lr_slab_{output_suffix}.log')}")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(os.path.join(slab_dir, f'lr_slab_{output_suffix}.tsv'), sep='\t', index=False)
    print(f"saved data in {os.path.join(slab_dir, f'lr_slab_{output_suffix}.tsv')}")

    # Plot parity with color by 'row' and marker by 'coord'
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}

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
                marker=coords[coords['coord'] == c]['marker'].values[0],
            )
            for _, row_data in subset.iterrows():
                row_features = pd.DataFrame([row_data[args.X].values], columns=args.X)
                y_pred_single = model.predict(row_features)[0]
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single))
    
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(bbox_to_anchor=(-0.10, 1), loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(slab_dir, f'lr_slab_{output_suffix}.png'))
    plt.close()

    # Save covariance and correlation matrices
    df_metrics = pd.concat([Y, X], axis=1)  # Place Y first
    cor = df_metrics.corr()

    plt.figure(figsize=(6, 5))
    cmap_forward = plt.get_cmap('coolwarm')
    cmap_reverse = plt.get_cmap('coolwarm_r')
    new_colors = np.vstack((
        cmap_reverse(np.linspace(0, 1, 256)),
        cmap_forward(np.linspace(0, 1, 256)),
    ))
    new_cmap = LinearSegmentedColormap.from_list('coolwarm_double', new_colors)
    ax = sns.heatmap(cor, annot=True, fmt='.2f', cmap=new_cmap, 
                annot_kws={"size": 7},
                cbar=False, vmin=-1, vmax=1)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(slab_dir, f'lr_slab_{output_suffix}_correlation.png'), 
                bbox_inches='tight', dpi=300, transparent=True)
    plt.close()

def calculate_cfse(coord, d_electrons, row):
    # Basic splitting patterns (energy differences)
    splitting_patterns = {
        # Octahedral structure
        'octahedral': {
            'orbitals': ['dxy', 'dxz', 'dyz', 'dx2y2', 'dz2'],
            'degeneracy': [3, 2],  # t2g(3), eg(2)
            'energies': [-0.4, 0.6]  # t2g: -0.4Δo, eg: +0.6Δo
        },
        # Tetrahedral structure
        'tetrahedral': {
            'orbitals': ['dx2y2', 'dz2', 'dxy', 'dxz', 'dyz'],
            'degeneracy': [2, 3],  # e(2), t2(3)
            'energies': [-0.6, 0.4]  # e: -0.6Δt, t2: +0.4Δt
        },
        # Square planar structure
        'square_planar': {
            'orbitals': ['dx2y2', 'dxy', 'dz2', 'dxz', 'dyz'],
            'degeneracy': [2, 1, 1, 1],
            'energies': [-0.35, -0.3, 0.1, 0.9]
        },
        # Pyramidal structure
        'pyramidal': {
            'orbitals': ['dx2y2', 'dxy', 'dz2', 'dxz', 'dyz'],
            'degeneracy': [2, 1, 1, 1],
            'energies': [-0.35, -0.3, 0.1, 0.9]
        }
    }
    
    # Coordination mapping
    coord_mapping = {
        'RS': 'octahedral',
        '+3': 'octahedral',
        '+4': 'octahedral',
        '+5': 'octahedral',
        '+6': 'octahedral',
        'ZB': 'tetrahedral',
        'WZ': 'tetrahedral',
        'TN': 'square_planar',
        'PD': 'square_planar',
        'NB': 'square_planar',
        'LT': 'pyramidal',
        'WW': 'tetrahedral',
        'ZZ': 'tetrahedral',
    }
    
    struct_type = coord_mapping[coord]
    pattern = splitting_patterns[struct_type]
    electron_config = get_electron_configuration(d_electrons, pattern['degeneracy'], pattern['energies'], row)
    
    cfse = 0.0
    for n_electrons, energy in zip(electron_config, pattern['energies']):
        cfse += n_electrons * energy
    return cfse

def get_electron_configuration(n_electrons, degeneracies, energies, row):
    config = [0] * len(energies)
    if row in ['4d', '5d']:
        # octaheral일 경우 직접 d 전자 배치 처리
        if len(degeneracies) == 2 and degeneracies[0] == 3 and degeneracies[1] == 2:  # octahedral (t2g, eg)
            if n_electrons <= 6:
                config = [n_electrons, 0]
            else:
                config = [6, n_electrons - 6]
        # tetrahedral일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 2 and degeneracies[0] == 2 and degeneracies[1] == 3:  # tetrahedral (e, t2)
            if n_electrons <= 4:
                config = [n_electrons, 0]
            else:
                config = [4, n_electrons - 4]
        # square planar일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 4 and degeneracies[0] == 2 and degeneracies[1] == 1 and degeneracies[2] == 1 and degeneracies[3] == 1:
            if n_electrons <= 4:
                config = [n_electrons, 0, 0, 0]
            elif n_electrons <= 6:
                config = [4, n_electrons - 4, 0, 0]
            elif n_electrons <= 8:
                config = [4, 2, n_electrons - 6, 0]
            else:
                config = [4, 2, 2, n_electrons - 8]
        # 일반적인 high spin 배치 로직 (square planar 등)
        else:
            remaining = n_electrons
            for i, deg in enumerate(degeneracies):
                if remaining >= 2 * deg:
                    config[i] = 2 * deg  # fully paired
                    remaining -= 2 * deg
                elif remaining > 0:
                    config[i] = remaining
                    remaining = 0
                else:
                    break
    # High spin case: magnetic moment가 0보다 큰 경우
    else:
        # octahedral일 경우 직접 d 전자 배치 처리
        if len(degeneracies) == 2 and degeneracies[0] == 3 and degeneracies[1] == 2:  # octahedral (t2g, eg)
            if n_electrons <= 3:
                config = [n_electrons, 0]
            elif n_electrons <= 5:
                config = [3, n_electrons - 3]
            elif n_electrons <= 8:
                config = [n_electrons - 2, 2]
            elif n_electrons <= 10:
                config = [6, n_electrons - 6]
        # tetrahedral일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 2 and degeneracies[0] == 2 and degeneracies[1] == 3:  # tetrahedral (e, t2)
            if n_electrons <= 2:
                config = [n_electrons, 0]
            elif n_electrons <= 5:
                config = [2, n_electrons - 2]
            elif n_electrons <= 7:
                config = [n_electrons - 3, 3]
            elif n_electrons <= 10:
                config = [4, n_electrons - 4]
        # square planar일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 4 and degeneracies[0] == 2 and degeneracies[1] == 1 and degeneracies[2] == 1 and degeneracies[3] == 1:
            if n_electrons <= 2:
                config = [n_electrons, 0, 0, 0]
            elif n_electrons <= 3:
                config = [2, 1, 0, 0]
            elif n_electrons <= 4:
                config = [2, 1, 1, 0]
            elif n_electrons <= 5:
                config = [2, 1, 1, 1]    
            elif n_electrons <= 7:
                config = [n_electrons - 4, 2, 1, 1]
            elif n_electrons <= 8:
                config = [4, 2, 1, 1]
            elif n_electrons <= 9:
                config = [4, 2, 2, 1]
            elif n_electrons <= 10:
                config = [4, 2, 2, 2]
        # 일반적인 high spin 배치 로직 (square planar 등)
        else:
            remaining = n_electrons
            total_orbitals = sum(degeneracies)
            single_filled = min(remaining, total_orbitals)
            electrons_left = single_filled
            for i, deg in enumerate(degeneracies):
                config[i] = min(deg, electrons_left)
                electrons_left -= config[i]
            remaining -= single_filled
            if remaining > 0:
                for i, deg in enumerate(degeneracies):
                    already_filled = config[i]
                    can_add = min(remaining, deg - already_filled)
                    config[i] += can_add
                    remaining -= can_add
                    if remaining == 0:
                        break
    return config


if __name__ == '__main__':
    main()