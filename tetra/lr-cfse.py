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

coords_data = [
    {'coord': 'WZ', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '1_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
    {'coord': 'ZB', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '2_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
    {'coord': 'TN', 'CN': 4, 'OS': 2, 'MN': 4, 'coord_dir': '3_SquarePlanar_TN', 'zorder': 3, 'marker': 'o', 'color': 'dodgerblue',},
    {'coord': 'PD', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '4_SquarePlanar_PD', 'zorder': 2, 'marker': 'o', 'color': 'deepskyblue',},
    {'coord': 'NB', 'CN': 4, 'OS': 2, 'MN': 6, 'coord_dir': '5_SquarePlanar_NB', 'zorder': 1, 'marker': 's', 'color': 'limegreen',},
    {'coord': 'RS', 'CN': 6, 'OS': 2, 'MN': 2, 'coord_dir': '6_Octahedral_RS',   'zorder': 6, 'marker': 'd', 'color': 'orchid',},
    {'coord': 'LT', 'CN': 4, 'OS': 2, 'MN': 2, 'coord_dir': '7_Pyramidal_LT',    'zorder': 0, 'marker': 'h', 'color': 'silver',},
]

coords = pd.DataFrame(coords_data).set_index('coord')
coords.index.name = None

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

indice = [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])]

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
    
    return selected_features

# CFSE 관련 함수들 추가
def calculate_cfse(coord, d_electrons, field_strength, mag, row):
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
        'NB': 'square_planar'
    }
    
    if coord not in coord_mapping:
        return {
            'base_cfse': 0.0,
            'ee_repulsion': 0.0,
            'jt_effect': 0.0,
            'exchange_stabilization': 0.0,
            'total_cfse': 0.0
        }
    
    struct_type = coord_mapping[coord]
    pattern = splitting_patterns[struct_type]
    
    electron_config = calculate_electron_configuration(
        d_electrons, 
        pattern['degeneracy'], 
        pattern['energies'],
        mag,
        row
    )
    
    # Calculate base CFSE
    base_cfse = 0.0
    for n_electrons, energy in zip(electron_config, pattern['energies']):
        base_cfse += n_electrons * energy
    
    # Calculate exchange repulsion
    exchange_stabilization = calculate_exchange_stabilization(d_electrons)
    
    # Calculate electron-electron repulsion
    ee_repulsion = calculate_electron_repulsion(electron_config)
    
    # Calculate Jahn-Teller effect
    jt_effect = calculate_jahn_teller_effect(electron_config, struct_type)
    
    # Calculate final CFSE
    total_cfse = (base_cfse - ee_repulsion + jt_effect + exchange_stabilization) * field_strength
    
    return {
        'base_cfse': base_cfse,
        'ee_repulsion': ee_repulsion,
        'jt_effect': jt_effect,
        'exchange_stabilization': exchange_stabilization,
        'total_cfse': total_cfse
    }

def calculate_electron_configuration(n_electrons, degeneracies, energies, mag, row):
    config = [0] * len(energies)
    
    # magnetic moment가 NaN이거나 0인 경우 low spin configuration 가정
    if row in ['fm', '4d', '5d']:
        # octahedral일 경우 직접 d 전자 배치 처리
        if len(degeneracies) == 2 and degeneracies[0] == 3 and degeneracies[1] == 2:  # octahedral (t2g, eg)
            # t2g와 eg 오비탈에 전자를 채우는 논리
            if n_electrons <= 6:
                config = [n_electrons, 0]
            else:
                config = [6, n_electrons - 6]
        
        # tetrahedral일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 2 and degeneracies[0] == 2 and degeneracies[1] == 3:  # tetrahedral (e, t2)
            # e와 t2 오비탈에 전자를 채우는 논리
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
            # t2g와 eg 오비탈에 전자를 채우는 논리
            if n_electrons <= 3:  # 1-3개 전자: 모두 t2g에 배치
                config = [n_electrons, 0]
            elif n_electrons <= 5:  # 4-5개 전자: t2g에 3개, 나머지는 eg에 배치
                config = [3, n_electrons - 3]
            elif n_electrons <= 8:  # 6-8개 전자: t2g에 전자가 추가됨
                config = [n_electrons - 2, 2]
            elif n_electrons <= 10:  # 9-10개 전자: eg에 전자가 추가됨
                config = [6, n_electrons - 6]
        
        # tetrahedral일 경우 직접 d 전자 배치 처리
        elif len(degeneracies) == 2 and degeneracies[0] == 2 and degeneracies[1] == 3:  # tetrahedral (e, t2)
            # e와 t2 오비탈에 전자를 채우는 논리
            if n_electrons <= 2:  # 1-2개 전자: 모두 e에 배치
                config = [n_electrons, 0]
            elif n_electrons <= 5:  # 3-5개 전자: e에 2개, 나머지는 t2에 배치
                config = [2, n_electrons - 2]
            elif n_electrons <= 7:  # 6-7개 전자: t2에 전자가 추가됨
                config = [n_electrons - 3, 3]
            elif n_electrons <= 10:  # 8-10개 전자: e에 전자가 추가됨
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
            # 먼저 모든 오비탈에 하나씩 채움
            remaining = n_electrons
            total_orbitals = sum(degeneracies)
            single_filled = min(remaining, total_orbitals)
            
            # 각 오비탈에 하나씩 전자 배치
            electrons_left = single_filled
            for i, deg in enumerate(degeneracies):
                config[i] = min(deg, electrons_left)
                electrons_left -= config[i]
            
            # 남은 전자 배치 (낮은 에너지 오비탈부터)
            remaining -= single_filled
            if remaining > 0:
                for i, deg in enumerate(degeneracies):
                    # 이미 배치된 전자 수
                    already_filled = config[i]
                    # 추가로 배치할 수 있는 전자 수
                    can_add = min(remaining, deg - already_filled)
                    config[i] += can_add
                    remaining -= can_add
                    
                    if remaining == 0:
                        break
    
    return config

def calculate_exchange_stabilization(d_electrons):
    """Calculate exchange repulsion considering paired electrons"""
    exchange_stabilization = 0.0
    if d_electrons <= 5 :
        exchange_stabilization = -d_electrons * (d_electrons - 1) / 2
    else:
        exchange_stabilization = -(10 - d_electrons) * (9 - d_electrons) / 2
    return exchange_stabilization

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
        'Ca': 0, 'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 6, 'Co': 7, 'Ni': 8, 'Cu': 9, 'Zn': 10, 'Ga': 11, 'Ge': 12, 
        'Sr': 0, 'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 6, 'Rh': 7, 'Pd': 8, 'Ag': 9, 'Cd': 10, 'In': 11, 'Sn': 12, 
        'Ba': 0, 'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 6, 'Ir': 7, 'Pt': 8, 'Au': 9, 'Hg': 10, 'Tl': 11, 'Pb': 12
    }
    
    if pd.notna(oxidation_state) and metal is not None:
        n_electrons = base_d_electrons.get(metal, 0) + 2 - oxidation_state
        d_electrons = max(0, min(10, round(n_electrons)))
        return n_electrons, d_electrons
    
    return np.nan, np.nan

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
    df[['n_electrons', 'd_electrons']] = df.apply(
        lambda row: pd.Series(estimate_d_electrons(
            row['chg'] if 'chg' in df.columns else None,
            row['OS'] if 'OS' in df.columns else None,
            row['metal'] if 'metal' in df.columns else None
        )), 
        axis=1
    )
    
    # Calculate field_strength
    df['field_strength'] = df.apply(
        lambda row: estimate_field_strength(
            row['madelung'] if 'madelung' in df.columns else None,
            row['l_bond'] if 'l_bond' in df.columns else None,
            row['-ICOHP'] if '-ICOHP' in df.columns else None
        ), 
        axis=1
    )
    
    # Calculate all CFSE components
    cfse_components = df.apply(
        lambda row: calculate_cfse(
            row['coord'],
            row['d_electrons'],
            row['field_strength'],
            row['mag'],  # magnetic moment 추가
            row['row']  # row 추가
        ),
        axis=1
    )
    
    # Add all CFSE components as separate features
    df['base_cfse'] = cfse_components.apply(lambda x: x['base_cfse'])
    df['ee_repulsion'] = cfse_components.apply(lambda x: x['ee_repulsion'])
    df['jt_effect'] = cfse_components.apply(lambda x: x['jt_effect'])
    df['exchange_stabilization'] = cfse_components.apply(lambda x: x['exchange_stabilization'])
    # df['cfse'] = cfse_components.apply(lambda x: x['total_cfse'])
    df['cfse'] = df['base_cfse']*df['OS']
    
    return df

def plot_cfse_by_metal_row(df, save_path):
    cfse_ylabels = {
        'd_electrons': 'Number of d-electrons',
        'base_cfse': 'Base CFSE (eV)',
        'ee_repulsion': 'Electron-Electron Repulsion (eV)',
        'jt_effect': 'Jahn-Teller Effect (eV)',
        'exchange_stabilization': 'Exchange Stabilization (eV)',
        'field_strength': 'Crystal Field Strength',
        'cfse': 'Total CFSE (eV)'
    }
    
    for row in ['fm', '3d', '4d', '5d']:
        for feature in cfse_ylabels.keys():
            plt.figure(figsize=(8, 6))
            # 각 배위구조별로 플롯
            for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']:
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                
                # coords에서 스타일 정보 가져오기
                zorder = coords.loc[coord, 'zorder']
                marker = coords.loc[coord, 'marker']
                color = coords.loc[coord, 'color']
                
                plt.plot(subset['numb'], subset[feature], 
                        marker=marker, color=color,
                        linestyle='-', label=coord, zorder=zorder)
            
            plt.xticks(np.arange(len(metals[row])), metals[row])
            plt.xlabel("Metal Index")
            plt.ylabel(cfse_ylabels[feature])
            plt.legend()
            plt.tight_layout()
            
            png_name = f"cfse_{row}_{feature}.png"
            plt.savefig(os.path.join(save_path, png_name))
            plt.close()
            print(f"Figure saved as {png_name}")

def plot_cfse_by_coordination(df, save_path):
    cfse_ylabels = {
        'd_electrons': 'Number of d-electrons',
        'base_cfse': 'Base CFSE (eV)',
        'ee_repulsion': 'Electron-Electron Repulsion (eV)',
        'jt_effect': 'Jahn-Teller Effect (eV)',
        'exchange_stabilization': 'Exchange Stabilization (eV)',
        'field_strength': 'Crystal Field Strength',
        'cfse': 'Total CFSE (eV)'
    }
    
    for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT']:
        for feature in cfse_ylabels.keys():
            plt.figure(figsize=(8, 6))
            
            # 스타일 정보 가져오기
            marker = coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
            colors = cmap(np.linspace(0.0, 0.6, 3))
            
            # 각 row별로 플롯
            for r, row in enumerate(['3d', '4d', '5d']):
                color = 'lightgray' if row == 'fm' else colors[r]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                plt.plot(subset['numb'], subset[feature],
                        marker=marker, color=color,
                        linestyle='-', label=row)
            
            plt.xticks(np.arange(len(indice)), indice)
            plt.xlabel("Metal Index")
            plt.ylabel(cfse_ylabels[feature])
            plt.legend()
            plt.tight_layout()
            
            png_name = f"cfse_{coord}_{feature}.png"
            plt.savefig(os.path.join(save_path, png_name))
            plt.close()
            print(f"Figure saved as {png_name}")

def group_highly_correlated_features(correlation_matrix, threshold=0.7):
    """
    상관계수의 절대값이 threshold보다 높은 피처들을 그룹화합니다.
    
    Args:
        correlation_matrix: 상관계수 행렬
        threshold: 그룹화 기준 상관계수 절대값 (default: 0.7)
    
    Returns:
        list: 그룹화된 피처들의 리스트
    """
    # 상관계수 행렬의 복사본 생성
    corr = correlation_matrix.copy()
    
    # 대각선 요소를 0으로 설정 (자기 자신과의 상관계수 제외)
    np.fill_diagonal(corr.values, 0)
    
    # 그룹화된 피처들을 저장할 리스트
    groups = []
    # 이미 그룹화된 피처들을 추적
    grouped_features = set()
    
    # 모든 피처에 대해 반복
    for feature in corr.columns:
        if feature in grouped_features:
            continue
            
        # 현재 피처와 상관계수가 threshold보다 높은 피처들 찾기
        high_corr = corr[abs(corr[feature]) > threshold].index.tolist()
        
        if high_corr:
            # 현재 피처를 포함한 그룹 생성
            group = [feature] + high_corr
            groups.append(group)
            # 그룹화된 피처들 기록
            grouped_features.update(group)
    
    # 그룹화되지 않은 단일 피처들을 각각 그룹으로 추가
    for feature in corr.columns:
        if feature not in grouped_features:
            groups.append([feature])
    
    return groups

def save_correlation_groups(groups, output_path):
    """
    그룹화된 피처들을 파일로 저장합니다.
    
    Args:
        groups: 그룹화된 피처들의 리스트
        output_path: 출력 파일 경로
    """
    with open(output_path, 'w') as f:
        f.write("Highly Correlated Feature Groups:\n\n")
        for i, group in enumerate(groups, 1):
            f.write(f"Group {i}:\n")
            f.write(", ".join(group) + "\n\n")

def main():
    parser = argparse.ArgumentParser(description='Linear regression using bulk_data.csv and mendeleev_data.csv')
    parser.add_argument('--Y', default='form', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=[
        'OS', 'CN', 'numb', 'group', 'mag', 'volume', 'l_bond', 'madelung',
        'chg', 'chgc', 'chgo', 'chgn',
        'ICOHP', 'ICOHPc', 'ICOHPo', 'ICOHPn',
        'ICOBI', 'ICOBIc', 'ICOBIo', 'ICOBIn',
        'ICOOP', 'ICOOPc', 'ICOOPo', 'ICOOPn',
        'ion-1', 'ion', 'ion+1', 'ionN-1', 'ionN', 'ionN+1', 
        'ion-1c', 'ionc', 'ion+1c', 'ionN-1c', 'ionNc', 'ionN+1c', 
        'ion-1o', 'iono', 'ion+1o', 'ionN-1o', 'ionNo', 'ionN+1o', 
        'ion-1n', 'ionn', 'ion+1n', 'ionN-1n', 'ionNn', 'ionN+1n', 
        'pauling', 'Natom', 'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
        'n_electrons', 'd_electrons', 'outer_e', 'base_cfse', 'ee_repulsion', 'jt_effect', 'field_strength', 'cfse', 'exchange_stabilization',
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

    # group 컬럼 추가
    df['group'] = df['numb'] + 3

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    # # Hf와 Ta의 +4 배위 구조 제외
    # df = df[~((df['metal'].isin(['Hf', 'Ta'])) & (df['coord'] == '+4'))]

    df['chgc'] = df['chg'] / df['CN']
    df['chgo'] = df['chg'] / df['OS']
    df['chgn'] = df['chg'] / df['OS'] / df['CN']
    
    # CFSE 피쳐 추가 - 위치 이동
    df = add_cfse_feature(df)
    
    # outer_e 계산 및 음수인 경우 제외
    df['outer_e'] = df['group'] - df['OS']
    df = df[df['outer_e'] >= 0]
    
    # Calculate ion values for each row
    df['ion-1'] = df.apply(lambda row: row[f'ion{int(row["OS"])-1}'], axis=1)
    df['ion'] = df.apply(lambda row: row[f'ion{int(row["OS"])}'], axis=1)
    df['ion+1'] = df.apply(lambda row: row[f'ion{int(row["OS"])+1}'], axis=1)

    df['ion-1o'] = df['ion-1'] / df['OS']
    df['iono'] = df['ion'] / df['OS']
    df['ion+1o'] = df['ion+1'] / df['OS']

    df['ion-1c'] = df['ion-1'] / df['CN']
    df['ionc'] = df['ion'] / df['CN']
    df['ion+1c'] = df['ion+1'] / df['CN']

    df['ion-1n'] = df['ion-1'] / df['OS'] / df['CN']
    df['ionn'] = df['ion'] / df['OS'] / df['CN']
    df['ion+1n'] = df['ion+1'] / df['OS'] / df['CN']

    df['-ICOHPo'] = df['-ICOHP'] / df['OS']
    df['ICOBIo'] = df['ICOBI'] / df['OS']
    df['-ICOOPo'] = df['-ICOOP'] / df['OS']

    df['-ICOHPc'] = df['-ICOHP'] / df['CN']
    df['ICOBIc'] = df['ICOBI'] / df['CN']
    df['-ICOOPc'] = df['-ICOOP'] / df['CN']

    df['-ICOHPn'] = df['-ICOHP'] / df['OS'] / df['CN']
    df['ICOBIn'] = df['ICOBI'] / df['OS'] / df['CN']
    df['-ICOOPn'] = df['-ICOOP'] / df['OS'] / df['CN']

    # Initialize columns with float type
    df['ionN-1'] = 0.0
    df['ionN'] = 0.0
    df['ionN+1'] = 0.0
    
    for idx, row in df.iterrows():
        for i in range(1, int(row['OS'])):
            df.at[idx, 'ionN-1'] = float(df.at[idx, 'ionN-1']) + float(row[f'ion{i}'])
        df.at[idx, 'ionN'] = float(df.at[idx, 'ionN-1']) + float(row[f'ion{int(row["OS"])}'])
        df.at[idx, 'ionN+1'] = float(df.at[idx, 'ionN']) + float(row[f'ion{int(row["OS"])+1}'])

    # ionN 관련 피쳐들을 OS로 나누어 정규화
    df['ionN-1o'] = df['ionN-1'] / df['OS']
    df['ionNo'] = df['ionN'] / df['OS']
    df['ionN+1o'] = df['ionN+1'] / df['OS']

    df['ionN-1c'] = df['ionN-1'] / df['CN']
    df['ionNc'] = df['ionN'] / df['CN']
    df['ionN+1c'] = df['ionN+1'] / df['CN']

    df['ionN-1n'] = df['ionN-1'] / df['OS'] / df['CN']
    df['ionNn'] = df['ionN'] / df['OS'] / df['CN']
    df['ionN+1n'] = df['ionN+1'] / df['OS'] / df['CN']

    # # 음수 버전의 피처 추가
    # negative_features = []
    # features_to_remove = []
    # for xx in ['Hevap', 'Hfus', 'Hform', 'ion-1', 'ion', 'ion+1', 'ion-1o', 'iono', 'ion+1o', 'ion-1c', 'ionc', 'ion+1c', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1', 'ionN-1o', 'ionNo', 'ionN+1o', 'ionN-1c', 'ionNc', 'ionN+1c', 'ionN-1n', 'ionNn', 'ionN+1n']:
    #     if xx in args.X:
    #         features_to_remove.append(xx)
    #         negative_features.append(f'-{xx}')
    
    # # 기존 피처 제거 및 음수 버전 추가
    # for feature in features_to_remove:
    #     args.X.remove(feature)
    # args.X.extend(negative_features)

    # # DataFrame 컬럼 수정
    # for col in ['Hevap', 'Hfus', 'Hform']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ion-1', 'ion', 'ion+1']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ion-1o', 'iono', 'ion+1o']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ion-1c', 'ionc', 'ion+1c']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ion-1n', 'ionn', 'ion+1n']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ionN-1', 'ionN', 'ionN+1']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ionN-1o', 'ionNo', 'ionN+1o']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ionN-1c', 'ionNc', 'ionN+1c']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # for col in ['ionN-1n', 'ionNn', 'ionN+1n']:
    #     if col in df.columns:
    #         df[f'-{col}'] = df[col] * -1

    # Drop rows with NaN in any relevant column
    df = df.dropna(subset=args.X + [args.Y])

    # Exclude ion1~7 columns from saving
    columns_to_save = [col for col in ['metal', 'row', 'coord'] + args.X + [args.Y] if col not in ['ion1', 'ion2', 'ion3', 'ion4', 'ion5', 'ion6', 'ion7']]
    df_to_save = df[columns_to_save]
    if args.output == 'result':
        df_to_save.to_csv(f'{root}/bulk_data_cfse.csv', sep=',')
        df_to_save.to_csv(f'{root}/bulk_data_cfse.tsv', sep='\t', float_format='%.2f')
    else:
        df_to_save.to_csv(f'{root}/bulk_data_cfse_{args.output}.csv', sep=',')
        df_to_save.to_csv(f'{root}/bulk_data_cfse_{args.output}.tsv', sep='\t', float_format='%.2f')

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
    with open(os.path.join(root, f'cfse_lr_{output_suffix}.log'), 'w') as f:
        f.write(f"Intercept: {model.intercept_:.4f}\n")
        for name, coef in zip(args.X, model.coef_):
            f.write(f"{name}: {coef:.4f}\n")
        f.write(f"\nR2: {r2:.4f}\nMAE: {mae:.4f}\nMSE: {mse:.4f}\n")

    # Save data with predictions
    df_result = df[['metal', 'row', 'coord'] + args.X].copy()
    df_result['Y_true'] = Y
    df_result['Y_pred'] = Y_pred
    df_result['residual'] = Y - Y_pred
    df_result.to_csv(os.path.join(root, f'cfse_lr_{output_suffix}.tsv'), sep='\t', index=False)

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
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single))

    
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(bbox_to_anchor=(-0.10, 1), loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(root, f'cfse_lr_{output_suffix}.png'))
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
    with open(os.path.join(root, f'cfse_lr_{args.output}_all.log'), 'w') as f:
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
    df_result_all.to_csv(os.path.join(root, f'cfse_lr_{args.output}_all.tsv'), sep='\t', index=False)

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
                plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single))

    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(bbox_to_anchor=(-0.10, 1), loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(root, f'cfse_lr_{args.output}_all.png'))
    plt.close()

    cmap_forward = plt.get_cmap('coolwarm')
    cmap_reverse = plt.get_cmap('coolwarm_r')

    new_colors = np.vstack((
        cmap_reverse(np.linspace(0, 1, 256)),
        cmap_forward(np.linspace(0, 1, 256)),
    ))

    new_cmap = LinearSegmentedColormap.from_list('coolwarm_double', new_colors)

    # Save covariance and correlation matrices
    cov.to_csv(os.path.join(root, f'cfse_covariance_{args.output}_all.tsv'), sep='\t')
    cor.to_csv(os.path.join(root, f'cfse_correlation_{args.output}_all.tsv'), sep='\t')

    # 다양한 threshold 값에 대해 그룹화 수행
    thresholds = [0.5, 0.6, 0.7, 0.8, 0.9]
    for threshold in thresholds:
        # 상관계수가 높은 피처들 그룹화
        correlation_groups = group_highly_correlated_features(cor, threshold)
        
        # 결과를 log 파일로 저장
        with open(os.path.join(root, f'cfse_correlation_groups_{args.output}_threshold_{threshold}.log'), 'w') as f:
            f.write(f"Correlation Threshold: {threshold}\n")
            f.write(f"Number of Groups: {len(correlation_groups)}\n\n")
            
            # 다중 피처 그룹과 단일 피처 그룹을 구분하여 출력
            multi_feature_groups = [g for g in correlation_groups if len(g) > 1]
            single_feature_groups = [g for g in correlation_groups if len(g) == 1]
            
            f.write("Multi-feature Groups:\n")
            for i, group in enumerate(multi_feature_groups, 1):
                f.write(f"Group {i}:\n")
                f.write(", ".join(group) + "\n\n")
            
            f.write("\nSingle-feature Groups:\n")
            single_features = [group[0] for group in single_feature_groups]
            f.write(", ".join(single_features) + "\n\n")
            
            f.write(f"\nTotal Groups: {len(correlation_groups)}\n")
            f.write(f"Multi-feature Groups: {len(multi_feature_groups)}\n")
            f.write(f"Single-feature Groups: {len(single_feature_groups)}\n")

    # Correlation matrix plot 수정
    plt.figure(figsize=(6, 5))
    ax = sns.heatmap(cor, annot=True, fmt='.2f', cmap=new_cmap, 
                annot_kws={"size": 7},
                cbar=False, vmin=-1, vmax=1)
    
    # x축 레이블 수정
    labels = ['Formation energy', 'Group number', 'Outer electron', 'Oxidation state', 'Coordination number', 'Ionization energy', 'Heat of evaporation', 'CFSE', '-ICOHP', 'Bader charge', 'Cell volume', 'M-O bond length', 'Magnetic moment']
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticklabels(labels, rotation=0)
    
    # 레이블 스타일 조정
    plt.tick_params(axis='both', which='major', labelsize=8)
    
    plt.tight_layout(pad=2.0)  # padding 증가
    plt.savefig(os.path.join(root, f'cfse_correlation_{args.output}_all.png'), 
                bbox_inches='tight', dpi=300, transparent=True)
    plt.close()

    # Covariance matrix plot도 동일하게 수정
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(cov, annot=True, fmt='.2f', cmap=new_cmap,
                cbar=False, vmin=-1, vmax=1)
    
    # x축 레이블 수정
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=12)
    ax.set_yticklabels(labels, rotation=0, fontsize=12)
    
    # 레이블 스타일 조정
    plt.tick_params(axis='both', which='major', labelsize=8)
    
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(root, f'cfse_covariance_{args.output}_all.png'),
                bbox_inches='tight', dpi=300)
    plt.close()

    print(f"Saved results using all features")

    # # Test different correlation thresholds
    # thresholds = [0.7, 0.8, 0.9]
    # threshold_suffixes = ['7', '8', '9']

    # for threshold, suffix in zip(thresholds, threshold_suffixes):
    #     print(f"\nTesting correlation threshold: {threshold}")
        
    #     # Select features based on correlation threshold
    #     selected_features = select_features_by_correlation(cor, args.Y, threshold)
        
    #     # Train new model with selected features
    #     X_selected = df[selected_features].astype(float)
    #     model_selected = LinearRegression()
    #     model_selected.fit(X_selected, Y)
        
    #     Y_pred_selected = model_selected.predict(X_selected)
    #     mae_selected = mean_absolute_error(Y, Y_pred_selected)
    #     mse_selected = mean_squared_error(Y, Y_pred_selected)
    #     r2_selected = model_selected.score(X_selected, Y)
        
    #     # Save results for selected features
    #     with open(os.path.join(root, f'cfse_lr_{args.output}_{suffix}.log'), 'w') as f:
    #         f.write(f"Correlation threshold: {threshold}\n")
    #         f.write(f"Selected Features: {', '.join(selected_features)}\n")
    #         f.write(f"Intercept: {model_selected.intercept_:.4f}\n")
    #         for name, coef in zip(selected_features, model_selected.coef_):
    #             f.write(f"{name}: {coef:.4f}\n")
    #         f.write(f"\nR2: {r2_selected:.4f}\nMAE: {mae_selected:.4f}\nMSE: {mse_selected:.4f}\n")

    #     # Save data with predictions
    #     df_result_selected = df[['metal', 'row', 'coord'] + selected_features].copy()
    #     df_result_selected['Y_true'] = Y
    #     df_result_selected['Y_pred'] = Y_pred_selected
    #     df_result_selected['residual'] = Y - Y_pred_selected
    #     df_result_selected.to_csv(os.path.join(root, f'cfse_lr_{args.output}_{suffix}.tsv'), sep='\t', index=False)

    #     plt.figure(figsize=(10, 8))
    #     for r in df['row'].unique():
    #         for c in df['coord'].unique():
    #             subset = df[(df['row'] == r) & (df['coord'] == c)]
    #             if subset.empty:
    #                 continue
    #             plt.scatter(
    #                 subset[args.Y],
    #                 model_selected.predict(subset[selected_features].astype(float)),
    #                 label=f'{r}_{c}',
    #                 alpha=0.3,
    #                 color=row_map.get(r, 'gray'),
    #                 marker=coord_map.get(c, 'x')
    #             )
    #             for _, row_data in subset.iterrows():
    #                 row_features = pd.DataFrame([row_data[selected_features].values], columns=selected_features)
    #                 y_pred_single = model_selected.predict(row_features)[0]
    #                 plt.annotate(row_data['metal'], (row_data[args.Y], y_pred_single))

    #     plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], '--', lw=1, color='black')
    #     plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    #     plt.ylabel(f'Predicted {ylabels[args.Y]}')
    #     plt.legend(loc='best')
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(root, f'cfse_lr_{args.output}_{suffix}.png'))
    #     plt.close()

    #     # Save covariance and correlation matrices for selected features
    #     df_metrics_selected = pd.concat([Y, X_selected], axis=1)  # Place Y first
    #     cov_selected = df_metrics_selected.cov()
    #     cor_selected = df_metrics_selected.corr()

    #     cov_selected.to_csv(os.path.join(root, f'cfse_covariance_{args.output}_{suffix}.tsv'), sep='\t')
    #     cor_selected.to_csv(os.path.join(root, f'cfse_correlation_{args.output}_{suffix}.tsv'), sep='\t')

    #     plt.figure(figsize=(18, 12))
    #     sns.heatmap(cov_selected, annot=True, fmt='.2f', cmap=new_cmap,
    #                 annot_kws={"size": 7},
    #                 cbar_kws={"shrink": 0.5, "aspect": 20, "pad": 0.02}, vmin=-1, vmax=1)
        
    #     plt.xticks(rotation=45, ha='right')
    #     plt.yticks(rotation=0)
        
    #     plt.tight_layout(pad=2.0)
    #     plt.savefig(os.path.join(root, f'cfse_covariance_{args.output}_{suffix}.png'),
    #                 bbox_inches='tight',
    #                 dpi=300)
    #     plt.close()

    #     plt.figure(figsize=(18, 12))
    #     sns.heatmap(cor_selected, annot=True, fmt='.2f', cmap=new_cmap,
    #                 annot_kws={"size": 7},
    #                 cbar_kws={"shrink": 0.5, "aspect": 20, "pad": 0.02}, vmin=-1, vmax=1)
        
    #     plt.xticks(rotation=45, ha='right')
    #     plt.yticks(rotation=0)
        
    #     plt.tight_layout(pad=2.0)
    #     plt.savefig(os.path.join(root, f'cfse_correlation_{args.output}_{suffix}.png'),
    #                 bbox_inches='tight',
    #                 dpi=300)
    #     plt.close()

    #     print(f"Saved results for threshold {threshold} with suffix {suffix}")

    # print(f"Saved: lr_{output_suffix}.log, lr_{output_suffix}.tsv")
    # print(f"Saved: covariance_{output_suffix}.tsv, covariance_{output_suffix}.png")
    # print(f"Saved: correlation_{output_suffix}.tsv, correlation_{output_suffix}.png")

    # CFSE 관련 플롯 생성
    # plot_cfse_by_metal_row(df, root)
    # plot_cfse_by_coordination(df, root)

if __name__ == '__main__':
    main()