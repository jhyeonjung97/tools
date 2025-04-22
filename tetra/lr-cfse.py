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

# CFSE 관련 함수들 추가
def calculate_cfse(coord, d_electrons, field_strength, mag):
    """Calculate Crystal Field Stabilization Energy with spin state consideration
    
    Args:
        coord (str): Coordination type
        d_electrons (int): Number of d electrons
        field_strength (float): Ligand field strength
        mag (float): Magnetic moment
    """
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
            'degeneracy': [1, 1, 1, 2],
            'energies': [0.9, 0.1, -0.3, -0.35]
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
        'NB': 'square_planar'
    }
    
    if coord not in coord_mapping:
        return {
            'base_cfse': 0.0,
            'ee_repulsion': 0.0,
            'jt_effect': 0.0,
            'total_cfse': 0.0
        }
        
    struct_type = coord_mapping[coord]
    pattern = splitting_patterns[struct_type]
    
    electron_config = calculate_electron_configuration(
        d_electrons, 
        pattern['degeneracy'], 
        pattern['energies'],
        mag
    )
    
    # Calculate base CFSE
    base_cfse = 0.0
    for n_electrons, energy in zip(electron_config, pattern['energies']):
        base_cfse += n_electrons * energy
    
    # Calculate electron-electron repulsion
    ee_repulsion = calculate_electron_repulsion(electron_config)
    
    # Calculate Jahn-Teller effect
    jt_effect = calculate_jahn_teller_effect(electron_config, struct_type)
    
    # Calculate final CFSE
    total_cfse = (base_cfse - ee_repulsion + jt_effect) * field_strength
    
    return {
        'base_cfse': base_cfse * field_strength,
        'ee_repulsion': ee_repulsion * field_strength,
        'jt_effect': jt_effect * field_strength,
        'total_cfse': total_cfse
    }

def calculate_electron_configuration(n_electrons, degeneracies, energies, mag):
    """Calculate electron configuration considering spin state based on magnetic moment
    
    Args:
        n_electrons (int): Total number of d electrons
        degeneracies (list): List of orbital degeneracies
        energies (list): List of orbital energies
        mag (float): Magnetic moment from DFT calculation
    
    Returns:
        list: Number of electrons in each orbital
    """
    config = [0] * len(energies)
    
    # magnetic moment가 NaN이거나 0인 경우 low spin configuration 가정
    if pd.isna(mag) or mag == 0:
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
        # magnetic moment로부터 unpaired electron 수 계산
        n_unpaired = round(float(mag))
        
        # 먼저 unpaired electron 배치
        remaining = n_electrons
        for i, deg in enumerate(degeneracies):
            if remaining >= deg:
                config[i] = min(deg, remaining)
                remaining -= deg
            else:
                config[i] = remaining
                break
    
    return config

def calculate_electron_repulsion(electron_config):
    """Calculate electron-electron repulsion considering paired electrons"""
    repulsion = 0.0
    for n_electrons in electron_config:
        if n_electrons > 1:
            # paired electron의 경우 반발이 더 큼
            n_pairs = n_electrons // 2
            repulsion += 0.2 * (n_electrons * (n_electrons - 1)) / 2
            repulsion += 0.1 * n_pairs  # additional repulsion for paired electrons
    return repulsion

def calculate_jahn_teller_effect(electron_config, struct_type):
    """Calculate Jahn-Teller effect"""
    if struct_type == 'octahedral':
        eg_electrons = electron_config[-2:]
        if eg_electrons[0] != eg_electrons[1]:
            return -0.1
    elif struct_type == 'tetrahedral':
        t2_electrons = electron_config[-3:]
        if any(n != t2_electrons[0] for n in t2_electrons):
            return -0.05
    
    return 0.0

def estimate_d_electrons(charge, oxidation_state, metal=None):
    """
    Estimate number of d electrons based on oxidation state
    
    Args:
        charge (float): Charge state from DFT calculation (not used)
        oxidation_state (float): Oxidation state (OS)
        metal (str): Metal symbol (e.g., 'Fe', 'Co', etc.)
        
    Returns:
        int: Estimated number of d electrons
    """
    base_d_electrons = {
        'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 6, 'Co': 7, 'Ni': 8, 'Cu': 9, 'Zn': 10,
        'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 6, 'Rh': 7, 'Pd': 8, 'Ag': 9, 'Cd': 10,
        'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 6, 'Ir': 7, 'Pt': 8, 'Au': 9, 'Hg': 10
    }
    
    if pd.notna(oxidation_state) and metal is not None:
        d_electrons = base_d_electrons.get(metal, 0) - oxidation_state
        return max(0, min(10, round(d_electrons)))
    
    return 0

def estimate_field_strength(madelung, bond_length, icohp):
    """Estimate ligand field strength"""
    if pd.isna(madelung) or pd.isna(bond_length) or pd.isna(icohp):
        return 0.0
        
    strength = (abs(madelung) * abs(icohp)) / (bond_length ** 5)
    return np.clip(strength, 0, 1)

def add_cfse_feature(df):
    """Add CFSE feature to dataframe considering magnetic moment"""
    df = df.copy()
    
    # Calculate d_electrons
    df['d_electrons'] = df.apply(
        lambda row: estimate_d_electrons(
            row['chg'] if 'chg' in df.columns else None,
            row['OS'] if 'OS' in df.columns else None,
            row['metal'] if 'metal' in df.columns else None
        ), 
        axis=1
    )
    
    # Calculate field_strength
    df['field_strength'] = df.apply(
        lambda row: estimate_field_strength(
            row['madelung'] if 'madelung' in df.columns else None,
            row['l_bond'] if 'l_bond' in df.columns else None,
            row['-ICOHPm'] if '-ICOHPm' in df.columns else None
        ), 
        axis=1
    )
    
    # Calculate all CFSE components
    cfse_components = df.apply(
        lambda row: calculate_cfse(
            row['coord'],
            row['d_electrons'],
            row['field_strength'],
            row['mag']  # magnetic moment 추가
        ),
        axis=1
    )
    
    # Add all CFSE components as separate features
    df['base_cfse'] = cfse_components.apply(lambda x: x['base_cfse'])
    df['ee_repulsion'] = cfse_components.apply(lambda x: x['ee_repulsion'])
    df['jt_effect'] = cfse_components.apply(lambda x: x['jt_effect'])
    df['cfse'] = cfse_components.apply(lambda x: x['total_cfse'])
    
    return df

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

    # CFSE 관련 피쳐들 추가
    cfse_features = ['base_cfse', 'ee_repulsion', 'jt_effect', 'field_strength', 'cfse']
    for feature in cfse_features:
        if feature not in args.X:
            args.X.append(feature)
    
    # CFSE 피쳐 추가 - 위치 이동
    df = add_cfse_feature(df)
    
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

    # Correlation matrix plot 수정
    plt.figure(figsize=(16, 12))  # 크기 증가
    sns.heatmap(cor, annot=True, fmt='.2f', cmap='coolwarm', 
                annot_kws={"size": 7},
                cbar_kws={"shrink": 0.5, "aspect": 20})
    
    # x축 레이블 수정
    plt.xticks(rotation=45, ha='right')  # rotation 각도 조정, 정렬 방식 변경
    # y축 레이블 수정
    plt.yticks(rotation=0)
    
    plt.tight_layout(pad=2.0)  # padding 증가
    plt.savefig(os.path.join(root, f'correlation_{args.output}_all.png'), 
                bbox_inches='tight',  # 여백 자동 조정
                dpi=300)  # 해상도 증가
    plt.close()

    # Covariance matrix plot도 동일하게 수정
    plt.figure(figsize=(16, 12))
    sns.heatmap(cov, annot=True, fmt='.2f', cmap='coolwarm',
                annot_kws={"size": 7},
                cbar_kws={"shrink": 0.5, "aspect": 20})
    
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(root, f'covariance_{args.output}_all.png'),
                bbox_inches='tight',
                dpi=300)
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

        plt.figure(figsize=(16, 12))
        sns.heatmap(cov_selected, annot=True, fmt='.2f', cmap='coolwarm',
                    annot_kws={"size": 7},
                    cbar_kws={"shrink": 0.5, "aspect": 20})
        
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(root, f'covariance_{args.output}_{suffix}.png'),
                    bbox_inches='tight',
                    dpi=300)
        plt.close()

        plt.figure(figsize=(16, 12))
        sns.heatmap(cor_selected, annot=True, fmt='.2f', cmap='coolwarm',
                    annot_kws={"size": 7},
                    cbar_kws={"shrink": 0.5, "aspect": 20})
        
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(root, f'correlation_{args.output}_{suffix}.png'),
                    bbox_inches='tight',
                    dpi=300)
        plt.close()

        print(f"Saved results for threshold {threshold} with suffix {suffix}")

    print(f"Saved: lr_{output_suffix}.log, lr_{output_suffix}.tsv")
    print(f"Saved: covariance_{output_suffix}.tsv, covariance_{output_suffix}.png")
    print(f"Saved: correlation_{output_suffix}.tsv, correlation_{output_suffix}.png")

if __name__ == '__main__':
    main()