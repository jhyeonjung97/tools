import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from tqdm import tqdm
import time
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic, WhiteKernel, ConstantKernel, RBF, Matern
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
from sklearn.model_selection import train_test_split, KFold, cross_val_score, GridSearchCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, max_error
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.utils.optimize import _check_optimize_result
import scipy.optimize
import socket
from sklearn.feature_selection import SequentialFeatureSelector, SelectFromModel
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.decomposition import PCA
import warnings
import scipy.stats as stats
warnings.filterwarnings('ignore')

# ANSI color codes
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'
CYAN = '\033[96m'
ENDC = '\033[0m'
BOLD = '\033[1m'

os.environ['JOBLIB_START_METHOD'] = 'fork'

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
    '-ICOHPm': '-ICOHP per Metal (eV)',
    'ICOBIm': 'ICOBI per Metal',
    '-ICOOPm': '-ICOOP per Metal (eV)',
    '-ICOHPn': '-ICOHP per Bond (eV)',
    'ICOBIn': 'ICOBI per Bond',
    '-ICOOPn': '-ICOOP per Bond (eV)',
    'madelung': 'Madelung Energy (Loewdin eV)',
}

def print_time(message, time_value):
    """Print time-related information with color coding"""
    print(f"{BOLD}{message}: {time_value:.2f} seconds{ENDC}")

class EnhancedGPR(GaussianProcessRegressor):
    def __init__(self, kernel=None, alpha=1e-10, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=0, 
                 normalize_y=False, copy_X_train=True, random_state=None, max_iter=1e06, gtol=1e-05):
        super().__init__(
            kernel=kernel,
            alpha=alpha,
            optimizer=optimizer,
            n_restarts_optimizer=n_restarts_optimizer,
            normalize_y=normalize_y,
            copy_X_train=copy_X_train,
            random_state=random_state
        )
        self.max_iter = max_iter
        self.gtol = gtol

    def _constrained_optimization(self, obj_func, initial_theta, bounds):
        if self.optimizer == "fmin_l_bfgs_b":
            opt_res = scipy.optimize.minimize(
                obj_func, 
                initial_theta, 
                method="L-BFGS-B", 
                jac=True, 
                bounds=bounds, 
                options={
                    'maxiter': self.max_iter,
                    'gtol': self.gtol,
                    'ftol': 1e-06,
                    'eps': 1e-08
                }
            )
            _check_optimize_result("lbfgs", opt_res)
            theta_opt, func_min = opt_res.x, opt_res.fun
        elif callable(self.optimizer):
            theta_opt, func_min = self.optimizer(obj_func, initial_theta, bounds=bounds)
        else:
            raise ValueError("Unknown optimizer %s." % self.optimizer)
        return theta_opt, func_min

def get_scaler(scaler_type):
    """Get appropriate scaler based on type"""
    if scaler_type == 'standard':
        return StandardScaler()
    elif scaler_type == 'robust':
        return RobustScaler()
    elif scaler_type == 'minmax':
        return MinMaxScaler()
    else:
        raise ValueError(f"Unknown scaler type: {scaler_type}")

def get_imputer(imputer_type):
    """Get appropriate imputer based on type"""
    if imputer_type == 'mean':
        return SimpleImputer(strategy='mean')
    elif imputer_type == 'median':
        return SimpleImputer(strategy='median')
    elif imputer_type == 'knn':
        return KNNImputer(n_neighbors=5)
    else:
        raise ValueError(f"Unknown imputer type: {imputer_type}")

def get_model(model_type, random_state=42, n_features=None):
    """Get the appropriate model based on the specified type"""
    if model_type == 'gpr':
        # 커널 파라미터 조정
        kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e3)) * RBF(
            length_scale=1.0,  # 단일 length_scale 사용
            length_scale_bounds=(1e-2, 1e2)  # 범위를 좁혀서 과적합 방지
        ) + WhiteKernel(
            noise_level=1.0,
            noise_level_bounds=(1e-2, 1e2)  # 노이즈 레벨 증가
        )
        
        return EnhancedGPR(
            kernel=kernel,
            random_state=random_state,
            n_restarts_optimizer=5,  # 재시작 횟수 감소
            alpha=1e-1,  # 정규화 강화
            optimizer='fmin_l_bfgs_b',
            normalize_y=True,
            max_iter=1000,  # 반복 횟수 제한
            gtol=1e-3  # 수렴 기준 완화
        )
    elif model_type == 'gbr':
        return GradientBoostingRegressor(
            n_estimators=300,
            learning_rate=0.01,
            max_depth=4,
            min_samples_split=10,
            min_samples_leaf=5,
            subsample=0.8,
            max_features=0.8,
            random_state=random_state
        )
    elif model_type == 'rf':
        return RandomForestRegressor(
            n_estimators=300,
            max_depth=None,
            min_samples_split=10,
            min_samples_leaf=5,
            max_features='auto',
            random_state=random_state
        )
    else:  # lr
        return Ridge(
            alpha=1.0,
            random_state=random_state
        )

def feature_selection(X, y, method='sequential', model_type='gbr', n_features=None):
    """Perform feature selection using different methods"""
    if method == 'sequential':
        selector = SequentialFeatureSelector(
            get_model(model_type),
            n_features_to_select=n_features,
            direction='forward',
            scoring='neg_mean_absolute_error',
            cv=5
        )
    elif method == 'model':
        selector = SelectFromModel(
            get_model(model_type),
            max_features=n_features
        )
    else:
        raise ValueError(f"Unknown feature selection method: {method}")
    
    selector.fit(X, y)
    return selector.get_support(indices=True)

def optimize_hyperparameters(model, X, y, model_type, cv=5):
    """Optimize hyperparameters using grid search"""
    if model_type == 'gpr':
        param_grid = {
            'kernel__k1__k2__length_scale': [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0],
            'kernel__k1__k2__alpha': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0],
            'kernel__k2__noise_level': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
        }
    elif model_type == 'gbr':
        param_grid = {
            'n_estimators': [100, 200, 300],
            'learning_rate': [0.01, 0.05, 0.1],
            'max_depth': [3, 4, 5],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4]
        }
    elif model_type == 'rf':
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20, 30],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4]
        }
    elif model_type == 'lr':
        param_grid = {
            'alpha': [0.1, 1.0, 10.0]
        }
    else:
        return model

    grid_search = GridSearchCV(
        model,
        param_grid,
        cv=cv,
        scoring='neg_mean_absolute_error',
        n_jobs=-1
    )
    grid_search.fit(X, y)
    return grid_search.best_estimator_

def get_prediction_std(model, X, model_type):
    """
    Calculate prediction standard deviation for different model types
    """
    if model_type == 'gpr':
        # For GPR, use the built-in uncertainty estimation
        _, std = model.predict(X, return_std=True)
        return std
    elif model_type == 'gbr':
        # For GBR, use the standard deviation of predictions from individual trees
        predictions = np.array([tree.predict(X) for tree in model.estimators_])
        return np.std(predictions, axis=0)
    elif model_type == 'rf':
        # For RF, use the standard deviation of predictions from individual trees
        predictions = np.array([tree.predict(X) for tree in model.estimators_])
        return np.std(predictions, axis=0)
    elif model_type == 'lr':
        # For LR, use the standard deviation of the residuals as a simple uncertainty estimate
        y_pred = model.predict(X)
        residuals = y_pred - model.predict(X)
        return np.std(residuals) * np.ones(len(X))
    else:
        # For other models, return zeros as a fallback
        return np.zeros(len(X))

def plot_feature_importance(X_train, X_test, y_train, y_test, feature_names, output_path, model_type):
    """Create a cumulative MAE plot for feature importance"""
    single_feature_scores = []
    cumulative_scores = []

    # Convert X_train and X_test to DataFrame if they are not already
    if isinstance(X_train, np.ndarray):
        X_train = pd.DataFrame(X_train, columns=feature_names)
    if isinstance(X_test, np.ndarray):
        X_test = pd.DataFrame(X_test, columns=feature_names)

    for feature in feature_names:
        # Select single feature using DataFrame indexing
        X_train_selected = X_train[[feature]].copy()
        X_test_selected = X_test[[feature]].copy()
        
        # Scale the selected features
        scaler_selected = RobustScaler()
        X_train_scaled_selected = scaler_selected.fit_transform(X_train_selected)
        X_test_scaled_selected = scaler_selected.transform(X_test_selected)
        
        # Train model with selected features
        if model_type == 'gpr':
            kernel = ConstantKernel(1.0) * RationalQuadratic() + WhiteKernel()
            model_selected = EnhancedGPR(kernel=kernel, random_state=42)
        else:  # gbr or lr
            model_selected = GradientBoostingRegressor(random_state=42)
        
        model_selected.fit(X_train_scaled_selected, y_train)
        
        # Make predictions
        if model_type == 'gpr':
            y_pred_train, _ = model_selected.predict(X_train_scaled_selected, return_std=True)
            y_pred_test, _ = model_selected.predict(X_test_scaled_selected, return_std=True)
        else:
            y_pred_train = model_selected.predict(X_train_scaled_selected)
            y_pred_test = model_selected.predict(X_test_scaled_selected)
        
        # Calculate MAE
        mae_train = mean_absolute_error(y_train, y_pred_train)
        mae_test = mean_absolute_error(y_test, y_pred_test)
        
        single_feature_scores.append({
            'feature': feature,
            'mae_train': mae_train,
            'mae_test': mae_test
        })
    
    # Sort features by individual MAE (best to worst - lower MAE means better performance)
    sorted_features = sorted(single_feature_scores, key=lambda x: x['mae_train'])
    sorted_feature_names = [score['feature'] for score in sorted_features]
    
    # Calculate cumulative MAE by adding one feature at a time (starting from most important)
    for i in range(1, len(sorted_feature_names) + 1):
        selected_features = sorted_feature_names[:i]
        X_train_selected = X_train[selected_features]
        X_test_selected = X_test[selected_features]
        
        # Scale the selected features
        scaler_selected = RobustScaler()
        X_train_scaled_selected = scaler_selected.fit_transform(X_train_selected)
        X_test_scaled_selected = scaler_selected.transform(X_test_selected)
        
        # Train model with selected features
        if model_type == 'gpr':
            kernel = ConstantKernel(1.0) * RationalQuadratic() + WhiteKernel()
            model_selected = EnhancedGPR(kernel=kernel, random_state=42)
        else:  # gbr or lr
            model_selected = GradientBoostingRegressor(random_state=42)
        
        model_selected.fit(X_train_scaled_selected, y_train)
        
        # Make predictions
        if model_type == 'gpr':
            y_pred_train, _ = model_selected.predict(X_train_scaled_selected, return_std=True)
            y_pred_test, _ = model_selected.predict(X_test_scaled_selected, return_std=True)
        else:
            y_pred_train = model_selected.predict(X_train_scaled_selected)
            y_pred_test = model_selected.predict(X_test_scaled_selected)
        
        # Calculate MAE
        mae_train = mean_absolute_error(y_train, y_pred_train)
        mae_test = mean_absolute_error(y_test, y_pred_test)
        
        cumulative_scores.append({
            'n_features': i,
            'features': selected_features,
            'mae_train': mae_train,
            'mae_test': mae_test
        })
    
    # Convert to DataFrame for plotting
    df_scores = pd.DataFrame(cumulative_scores)
    
    # Plot cumulative MAE vs number of features
    plt.figure(figsize=(10, 4))
    plt.plot(df_scores['n_features'], df_scores['mae_test'], '--o', color='blue', label='test')
    plt.plot(df_scores['n_features'], df_scores['mae_train'], '--o', color='orange', label='train')
    
    # Customize plot
    plt.xlabel('Number of descriptors')
    plt.ylabel('Formation Energy MAE (eV)')
    plt.legend()
    
    # Set x-axis limits with margins and minor ticks
    x_margin = 0.5
    plt.xlim(1 - x_margin, len(sorted_feature_names) + x_margin)
    ax = plt.gca()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    
    # Add feature names on top with matching limits
    ax2 = plt.gca().twiny()
    ax2.set_xlim(1 - x_margin, len(sorted_feature_names) + x_margin)
    ax2.set_xticks(range(1, len(sorted_feature_names) + 1))
    ax2.set_xticklabels(sorted_feature_names, rotation=90)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path)
    plt.close()
    print(f"{BLUE}Feature importance plot saved as {output_path}{ENDC}")

    return cumulative_scores

def get_preferred_coordination(energies):
    """Get the coordination with minimum energy"""
    min_energy = min(energies.values())
    return [coord for coord, energy in energies.items() 
            if np.isclose(energy, min_energy, rtol=1e-10, atol=1e-10)][0]

def get_coordination_type(coord):
    """Get the coordination type based on the coordination"""
    if coord in ['WZ', 'ZB']:
        return 'tetrahedral'
    elif coord in ['TN', 'PD', 'NB']:
        return 'squareplanar'
    elif coord == 'RS':
        return 'octahedral'
    elif coord == 'LT':
        return 'pyramidal'
    return 'unknown'

def analyze_coordination_preference(df_bulk, df_pred, energy_threshold=0.2, target='form'):
    """
    Compare DFT and predicted coordination preferences
    
    Args:
        df_bulk: DataFrame with DFT results
        df_pred: DataFrame with predictions
        energy_threshold: Energy threshold for considering multiple coordinations
        target: Target energy column ('form' or 'coh')
    
    Returns:
        DataFrame with comparison results
    """
    results = []
    
    # Group prediction data by metal and coordination
    pred_grouped = df_pred.groupby(['metal', 'coord'])['Y_pred'].mean().reset_index()
    pred_grouped = pred_grouped.pivot(index='metal', columns='coord', values='Y_pred')
    
    # Print only essential information
    bulk_metals = set(df_bulk['metal'].unique())
    pred_metals = set(pred_grouped.index)
    missing_metals = bulk_metals - pred_metals
    if missing_metals:
        print(f"{YELLOW}Metals without predictions: {missing_metals}{ENDC}")
    
    for metal in df_bulk['metal'].unique():
        try:
            # DFT data
            dft_data = df_bulk[df_bulk['metal'] == metal]
            if len(dft_data) == 0:
                continue
                
            dft_energies = dict(zip(dft_data['coord'], dft_data[target]))
            dft_energies = {k: v for k, v in dft_energies.items() if pd.notna(v)}  # Remove NaN values
            
            if not dft_energies:
                continue
                
            dft_min_energy = min(dft_energies.values())
            dft_preferred = get_preferred_coordination(dft_energies)
            dft_type = get_coordination_type(dft_preferred)
            
            # Predicted data
            if metal not in pred_grouped.index:
                continue
                
            pred_energies = pred_grouped.loc[metal].dropna().to_dict()
            if not pred_energies:
                continue
                
            pred_min_energy = min(pred_energies.values())
            pred_preferred = get_preferred_coordination(pred_energies)
            pred_type = get_coordination_type(pred_preferred)
            
            # Find all coordinations within threshold
            dft_preferred_all = [coord for coord, energy in dft_energies.items() 
                               if energy <= dft_min_energy + energy_threshold]
            pred_preferred_all = [coord for coord, energy in pred_energies.items() 
                                if energy <= pred_min_energy + energy_threshold]
            
            results.append({
                'metal': metal,
                'row': dft_data['row'].iloc[0],
                'dft_preferred': dft_preferred,
                'pred_preferred': pred_preferred,
                'dft_type': dft_type,
                'pred_type': pred_type,
                'match': dft_preferred == pred_preferred,
                'type_match': dft_type == pred_type,
                'dft_all_preferred': ', '.join(dft_preferred_all),
                'pred_all_preferred': ', '.join(pred_preferred_all),
                'dft_energies': dft_energies,
                'pred_energies': pred_energies,
                'energy_diff': abs(dft_min_energy - pred_min_energy)
            })
        except Exception as e:
            print(f"{RED}Error processing metal {metal}: {str(e)}{ENDC}")
            continue
    
    return pd.DataFrame(results)

def plot_coordination_comparison(results, output_dir, output_suffix, model_type, target, row_str):
    """Create visualization plots for the coordination comparison"""
    # 1. Heatmap showing matches and mismatches for specific coordinations
    plt.figure(figsize=(4, 3))
    
    # Create comparison matrix
    all_coords = ['ZB', 'WZ', 'RS', 'NB', 'PD', 'TN', 'LT']
    comparison_matrix = np.zeros((len(all_coords), len(all_coords)))
    coord_to_idx = {coord: i for i, coord in enumerate(all_coords)}
    
    for _, row in results.iterrows():
        dft_idx = coord_to_idx[row['dft_preferred']]
        pred_idx = coord_to_idx[row['pred_preferred']]
        comparison_matrix[dft_idx, pred_idx] += 1
    
    # Plot heatmap
    sns.heatmap(comparison_matrix, 
                xticklabels=all_coords, 
                yticklabels=all_coords,
                annot=True, 
                fmt='g',
                cmap='Reds',
                vmin=0)
    plt.xlabel('Predicted Preferred Coordination')
    plt.ylabel('DFT Preferred Coordination')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{target}_pred_{model_type}_{row_str}_{output_suffix}_coord.png'))
    plt.close()

    # 2. Heatmap showing matches and mismatches for coordination types
    plt.figure(figsize=(4, 3))
    
    # Create comparison matrix for types
    all_types = ['tetrahedral', 'squareplanar', 'octahedral', 'pyramidal']
    type_matrix = np.zeros((len(all_types), len(all_types)))
    type_to_idx = {coord_type: i for i, coord_type in enumerate(all_types)}
    
    for _, row in results.iterrows():
        dft_type_idx = type_to_idx[row['dft_type']]
        pred_type_idx = type_to_idx[row['pred_type']]
        type_matrix[dft_type_idx, pred_type_idx] += 1
    
    # Plot heatmap for types
    sns.heatmap(type_matrix, 
                xticklabels=all_types, 
                yticklabels=all_types,
                annot=True, 
                fmt='g',
                cmap='Reds',
                vmin=0)
    plt.xlabel('Predicted Coordination Type')
    plt.ylabel('DFT Coordination Type')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{target}_pred_{model_type}_{row_str}_{output_suffix}_type.png'))
    plt.close()

def remove_outliers(X, y, n_sigma=3):
    """Remove outliers based on z-score"""
    z_scores = np.abs(stats.zscore(y))
    mask = z_scores < n_sigma
    return X[mask], y[mask]

def main():
    start_time = time.time()
    print("Starting enhanced bulk prediction analysis...")
    
    parser = argparse.ArgumentParser(description='Enhanced bulk prediction with multiple optimization techniques')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr'], default='gbr',
                      help='Model type to use')
    parser.add_argument('--Y', type=str, choices=['form', 'coh'], default='form',
                      help='Target column from bulk_data.csv (form or coh)')
    parser.add_argument('--X', nargs='+', default=[
        'numb', 'chg', 'mag', 'volume', 'l_bond', 'n_bond',
        'grosspop', 'madelung', 'ICOHPm', 'ICOHPn', 'ICOBIm', 'ICOBIn', 'ICOOPm', 'ICOOPn', 
        'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
        'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size')
    parser.add_argument('--random_state', type=int, default=42, help='Random state')
    parser.add_argument('--scaler', type=str, choices=['standard', 'robust', 'minmax'], default='standard',
                      help='Type of scaler to use')
    parser.add_argument('--imputer', type=str, choices=['mean', 'median', 'knn'], default='mean',
                      help='Type of imputer to use')
    parser.add_argument('--feature_selection', type=str, choices=['sequential', 'model', 'none'], default='none',
                      help='Feature selection method')
    parser.add_argument('--n_features', type=int, default=None, help='Number of features to select')
    parser.add_argument('--pca', action='store_true', help='Use PCA for dimensionality reduction')
    parser.add_argument('--n_components', type=float, default=0.95, help='Variance ratio for PCA')
    parser.add_argument('--optimize', action='store_true', help='Optimize hyperparameters')
    parser.add_argument('--cv', type=int, default=5, help='Number of cross-validation folds')
    parser.add_argument('--threshold', type=float, default=0.2, help='Energy threshold for coordination analysis')
    args = parser.parse_args()
    
    # Convert feature names if they start with ICOHP or ICOOP (prepend '-')
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    
    # Get root directory
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
    
    # Define output paths
    output_suffix = args.output
    row_str = ''.join(sorted(args.row)) if args.row else 'all'
    log_path = os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}.log')
    tsv_path = os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}.tsv')
    png_path = os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}.png')
    importance_png_path = os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}_importance.png')
    json_path = os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}.json')

    try:
        # Load data
        print("Loading data...")
        load_start = time.time()
        df_bulk = pd.read_csv(os.path.join(root, 'bulk_data.csv'), index_col=0)
        df_mend = pd.read_csv(os.path.join(root, 'mendeleev_data.csv'), index_col=0)
        print_time("Data loading completed", time.time() - load_start)
    except FileNotFoundError as e:
        print(f"{RED}Error: Required data file not found: {e}{ENDC}")
        exit(1)

    print("Preprocessing data...")
    preprocess_start = time.time()
    
    # Merge bulk and mendeleev data
    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    df = df.rename(columns={'row_bulk': 'row', 'numb_mend': 'numb'})
    df = df.drop(columns=['row_mend', 'numb_bulk'])
    df = df[df['row'] != 'fm']

    # Handle missing values in features quietly
    for col in args.X + [args.Y]:  # Y값도 포함하여 NaN 처리
        if df[col].isna().any():
            if df[col].dtype in [np.float64, np.int64]:
                df[col] = df[col].fillna(df[col].median())
            else:
                df[col] = df[col].fillna(df[col].mode()[0])

    # Drop rows where y is NaN
    valid_indices = ~df[args.Y].isna()
    df = df[valid_indices]

    # Print dataset size
    print(f"\nDataset size before filtering:")
    print(f"Total rows: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    print(f"Features used: {len(args.X)}")
    print(f"Target variable: {args.Y}")

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    # Print dataset size after filtering
    print(f"\nDataset size after filtering:")
    print(f"Total rows: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    print(f"Features used: {len(args.X)}")
    print(f"Target variable: {args.Y}")

    # Handle missing values and data preparation
    print("Preprocessing data...")
    preprocess_start = time.time()

    # Handle missing values in features quietly
    for col in args.X + [args.Y]:  # Y값도 포함하여 NaN 처리
        if df[col].isna().any():
            if df[col].dtype in [np.float64, np.int64]:
                df[col] = df[col].fillna(df[col].median())
            else:
                df[col] = df[col].fillna(df[col].mode()[0])

    # Drop rows where y is NaN
    valid_indices = ~df[args.Y].isna()
    df = df[valid_indices]

    # Scale features
    print("Scaling features...")
    scale_start = time.time()
    scaler = RobustScaler()  # 또는 args.scaler에 따라 선택
    X = df[args.X].astype(float)
    y = df[args.Y].astype(float)
    
    # Remove outliers
    X_cleaned, y_cleaned = remove_outliers(X, y)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X_cleaned, y_cleaned, test_size=args.test_size, random_state=args.random_state
    )
    
    # Scale the features
    X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns, index=X_train.index)
    X_test_scaled = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns, index=X_test.index)
    print_time("Feature scaling completed", time.time() - scale_start)

    # Create df_train and df_test for metadata
    df_train = pd.DataFrame(X_train, columns=args.X)
    df_train['Y'] = y_train
    df_train['row'] = df.loc[X_train.index, 'row']
    df_train['coord'] = df.loc[X_train.index, 'coord']
    df_train['metal'] = df.loc[X_train.index, 'metal']
    
    df_test = pd.DataFrame(X_test, columns=args.X)
    df_test['Y'] = y_test
    df_test['row'] = df.loc[X_test.index, 'row']
    df_test['coord'] = df.loc[X_test.index, 'coord']
    df_test['metal'] = df.loc[X_test.index, 'metal']

    # Print training and test dataset sizes
    print(f"\nTraining and test dataset sizes:")
    print(f"Training data: {len(X_train)} samples")
    print(f"Test data: {len(X_test)} samples")
    print(f"Total data: {len(X_train) + len(X_test)} samples")

    # Feature selection
    if args.feature_selection != 'none':
        print("Performing feature selection...")
        feature_start = time.time()
        selected_indices = feature_selection(
            X_train_scaled, y_train, 
            method=args.feature_selection,
            model_type=args.model,
            n_features=args.n_features
        )
        X_train_scaled = X_train_scaled[:, selected_indices]
        X_test_scaled = X_test_scaled[:, selected_indices]
        selected_features = [args.X[i] for i in selected_indices]
        print(f"Selected features: {selected_features}")
        print_time("Feature selection completed", time.time() - feature_start)

    # PCA
    if args.pca:
        print("Applying PCA...")
        pca_start = time.time()
        pca = PCA(n_components=args.n_components)
        X_train_scaled = pca.fit_transform(X_train_scaled)
        X_test_scaled = pca.transform(X_test_scaled)
        print(f"Number of components after PCA: {X_train_scaled.shape[1]}")
        print_time("PCA completed", time.time() - pca_start)

    # Get and optimize model
    print(f"Setting up {args.model.upper()} model...")
    model = get_model(args.model, args.random_state)
    
    # 모델이 GPR인 경우 커널 파라미터 업데이트
    if args.model == 'gpr':
        # 각 특성별 길이 스케일로 커널 업데이트
        kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e3)) * RBF(
            length_scale=[1.0] * X_train.shape[1],
            length_scale_bounds=(1e-2, 1e2)
        ) + WhiteKernel(
            noise_level=1.0,
            noise_level_bounds=(1e-2, 1e2)
        )
        model.kernel = kernel
    
    # 교차 검증으로 과적합 모니터링
    cv = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=cv, scoring='neg_mean_squared_error')
    print(f"Cross-validation MSE: {-cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")

    # 모델 학습
    model.fit(X_train_scaled, y_train)

    # Make predictions
    print("Making predictions...")
    pred_start = time.time()
    y_pred_train = model.predict(X_train_scaled)
    y_pred_test = model.predict(X_test_scaled)
    print_time("Predictions completed", time.time() - pred_start)

    # Calculate metrics
    print("Calculating metrics...")
    metrics_start = time.time()
    metrics = {
        'train': {
            'r2': r2_score(y_train, y_pred_train),
            'mae': mean_absolute_error(y_train, y_pred_train),
            'mse': mean_squared_error(y_train, y_pred_train),
            'rmse': np.sqrt(mean_squared_error(y_train, y_pred_train)),
            'max_error': max_error(y_train, y_pred_train)
        },
        'test': {
            'r2': r2_score(y_test, y_pred_test),
            'mae': mean_absolute_error(y_test, y_pred_test),
            'mse': mean_squared_error(y_test, y_pred_test),
            'rmse': np.sqrt(mean_squared_error(y_test, y_pred_test)),
            'max_error': max_error(y_test, y_pred_test)
        }
    }
    print_time("Metrics calculation completed", time.time() - metrics_start)

    # Save results
    print("Saving results...")
    save_start = time.time()

    # Save metrics
    with open(log_path, 'w') as f:
        f.write("Training Metrics:\n")
        f.write(f"R2: {metrics['train']['r2']:.4f}\n")
        f.write(f"MAE: {metrics['train']['mae']:.4f}\n")
        f.write(f"MSE: {metrics['train']['mse']:.4f}\n")
        f.write(f"RMSE: {metrics['train']['rmse']:.4f}\n")
        f.write(f"Max Error: {metrics['train']['max_error']:.4f}\n\n")
        f.write("Test Metrics:\n")
        f.write(f"R2: {metrics['test']['r2']:.4f}\n")
        f.write(f"MAE: {metrics['test']['mae']:.4f}\n")
        f.write(f"MSE: {metrics['test']['mse']:.4f}\n")
        f.write(f"RMSE: {metrics['test']['rmse']:.4f}\n")
        f.write(f"Max Error: {metrics['test']['max_error']:.4f}\n")

    # Save predictions
    df_result = pd.DataFrame({
        'metal': df.loc[valid_indices, 'metal'],
        'row': df.loc[valid_indices, 'row'],
        'coord': df.loc[valid_indices, 'coord'],
        'Y_true': y,
        'Y_pred': model.predict(scaler.transform(X_cleaned))
    })
    df_result.to_csv(tsv_path, sep='\t', index=False)

    # Plot parity plot with uncertainty
    print("Creating parity plot...")
    plt.figure(figsize=(10, 8))

    # Create DataFrames for train and test sets
    train_indices = np.arange(len(X_train))
    test_indices = np.arange(len(X_test))
    
    df_train = pd.DataFrame(X_train_scaled, columns=args.X)
    df_train['Y'] = y_train.values
    df_train['Y_pred'] = y_pred_train
    df_train['std'] = get_prediction_std(model, X_train_scaled, model_type=args.model)
    df_train['row'] = df.iloc[train_indices]['row'].values
    df_train['coord'] = df.iloc[train_indices]['coord'].values
    df_train['metal'] = df.iloc[train_indices]['metal'].values

    df_test = pd.DataFrame(X_test_scaled, columns=args.X)
    df_test['Y'] = y_test.values
    df_test['Y_pred'] = y_pred_test
    df_test['std'] = get_prediction_std(model, X_test_scaled, model_type=args.model)
    df_test['row'] = df.iloc[test_indices]['row'].values
    df_test['coord'] = df.iloc[test_indices]['coord'].values
    df_test['metal'] = df.iloc[test_indices]['metal'].values

    # Define color and marker mappings
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}
    coord_map = {'WZ': '>', 'ZB': '<', 'TN': 'o', 'PD': 'o', 'NB': 's', 'RS': 'd', 'LT': 'h'}

    # Plot training data with error bars
    for r in df.iloc[train_indices]['row'].unique():
        for c in df.iloc[train_indices]['coord'].unique():
            # Plot training data
            subset_train = df_train[(df_train['row'] == r) & (df_train['coord'] == c)]
            if not subset_train.empty:
                plt.errorbar(
                    subset_train['Y'].values,
                    subset_train['Y_pred'].values,
                    yerr=subset_train['std'].values,
                    fmt=coord_map.get(c, 'x'),
                    label=f'{r}_{c}',
                    alpha=0.3,
                    color=row_map.get(r, 'gray'),
                    markeredgecolor=row_map.get(r, 'gray'),
                    markerfacecolor=row_map.get(r, 'gray'),
                    ecolor='silver',
                    capsize=0,
                    linewidth=0.5
                )
                # Add labels for training set
                for i, (_, row_data) in enumerate(subset_train.iterrows()):
                    plt.annotate(row_data['metal'], 
                               (row_data['Y'], row_data['Y_pred']),
                               fontsize=8)
            
            # Plot test data
            subset_test = df_test[(df_test['row'] == r) & (df_test['coord'] == c)]
            if not subset_test.empty:
                plt.errorbar(
                    subset_test['Y'].values,
                    subset_test['Y_pred'].values,
                    yerr=subset_test['std'].values,
                    fmt=coord_map.get(c, 'x'),
                    alpha=0.3,
                    color=row_map.get(r, 'gray'),
                    markeredgecolor=row_map.get(r, 'gray'),
                    markerfacecolor=row_map.get(r, 'gray'),
                    ecolor='silver',
                    capsize=0,
                    linewidth=0.5
                )
                # Add labels for test set
                for i, (_, row_data) in enumerate(subset_test.iterrows()):
                    plt.annotate(row_data['metal'], 
                               (row_data['Y'], row_data['Y_pred']),
                               fontsize=8)

    plt.plot([y.min(), y.max()], [y.min(), y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    
    # Save plot
    plt.savefig(png_path)
    plt.close()
    print(f"{BLUE}Parity plot saved as {png_path}{ENDC}")

    # Calculate feature importance
    print("Calculating feature importance...")
    importance_start = time.time()
    importance_scores = plot_feature_importance(
        X_train_scaled, X_test_scaled, y_train, y_test, 
        X_train_scaled.columns.tolist(), 
        os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}_importance.png'),
        args.model
    )
    print_time("Feature importance calculation completed", time.time() - importance_start)
    
    # Update results dictionary with feature importance
    results = {
        'metrics': metrics,
        'model_params': str(model.get_params()) if hasattr(model, 'get_params') else None,
        'feature_names': args.X,
        'target': args.Y,
        'preprocessing': {
            'scaler': args.scaler,
            'imputer': args.imputer,
            'feature_selection': args.feature_selection,
            'pca': args.pca,
            'n_components': args.n_components if args.pca else None
        },
        'feature_importance': importance_scores
    }
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)

    # Analyze coordination preferences
    print("Analyzing coordination preferences...")
    coord_results = analyze_coordination_preference(df, df_result, args.threshold, args.Y)
    
    if len(coord_results) > 0:
        # Save coordination comparison results
        coord_results.to_csv(os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}_coord.csv'))
        
        # Create coordination summary file
        with open(os.path.join(root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}_coord.log'), 'w') as f:
            # Overall statistics
            f.write("Overall Statistics:\n")
            f.write("-----------------\n")
            f.write(f"Total metals analyzed: {len(coord_results)}\n")
            if 'match' in coord_results.columns:
                correct_predictions = coord_results['match'].sum()
                incorrect_predictions = len(coord_results) - correct_predictions
                f.write(f"Correct predictions: {int(correct_predictions)} metals\n")
                f.write(f"Incorrect predictions: {int(incorrect_predictions)} metals\n")
                f.write(f"Overall accuracy: {coord_results['match'].mean():.2f}\n")
                f.write(f"Average energy difference: {coord_results['energy_diff'].mean():.2f} eV\n\n")
                
                # Add coordination type statistics
                type_correct = coord_results['type_match'].sum()
                type_incorrect = len(coord_results) - type_correct
                f.write("\nCoordination Type Statistics:\n")
                f.write("-------------------------\n")
                f.write(f"Correct type predictions: {int(type_correct)} metals\n")
                f.write(f"Incorrect type predictions: {int(type_incorrect)} metals\n")
                f.write(f"Type prediction accuracy: {coord_results['type_match'].mean():.2f}\n\n")
                
                # Add statistics by coordination type
                f.write("Statistics by Coordination Type:\n")
                f.write("-----------------------------\n")
                for coord_type in ['tetrahedral', 'squareplanar', 'octahedral', 'pyramidal']:
                    type_data = coord_results[coord_results['dft_type'] == coord_type]
                    if len(type_data) > 0:
                        f.write(f"\n{coord_type.capitalize()} coordination:\n")
                        f.write(f"Number of metals: {len(type_data)}\n")
                        type_correct = type_data['type_match'].sum()
                        type_incorrect = len(type_data) - type_correct
                        f.write(f"Correct type predictions: {int(type_correct)} metals\n")
                        f.write(f"Incorrect type predictions: {int(type_incorrect)} metals\n")
                        f.write(f"Type prediction accuracy: {type_data['type_match'].mean():.2f}\n")
                        f.write(f"Average energy difference: {type_data['energy_diff'].mean():.2f} eV\n")
            
            # Statistics by row
            f.write("\nStatistics by Row:\n")
            f.write("----------------\n")
            for row in ['3d', '4d', '5d']:
                row_data = coord_results[coord_results['row'] == row]
                f.write(f"\n{row} metals:\n")
                f.write(f"Number of metals: {len(row_data)}\n")
                if 'match' in coord_results.columns:
                    row_correct = row_data['match'].sum()
                    row_incorrect = len(row_data) - row_correct
                    f.write(f"Correct predictions: {int(row_correct)} metals\n")
                    f.write(f"Incorrect predictions: {int(row_incorrect)} metals\n")
                    f.write(f"Accuracy: {row_data['match'].mean():.2f}\n")
                    f.write(f"Average energy difference: {row_data['energy_diff'].mean():.2f} eV\n")
                    
                    # Add row-wise coordination type statistics
                    row_type_correct = row_data['type_match'].sum()
                    row_type_incorrect = len(row_data) - row_type_correct
                    f.write(f"Correct type predictions: {int(row_type_correct)} metals\n")
                    f.write(f"Incorrect type predictions: {int(row_type_incorrect)} metals\n")
                    f.write(f"Type prediction accuracy: {row_data['type_match'].mean():.2f}\n")
            
            # Detailed analysis of mismatches
            if 'match' in coord_results.columns:
                f.write("\nMismatch Analysis:\n")
                f.write("----------------\n")
                mismatches = coord_results[~coord_results['match']]
                for _, row in mismatches.iterrows():
                    f.write(f"\n{row['metal']} ({row['row']}):\n")
                    f.write(f"DFT preferred: {row['dft_preferred']} ({row['dft_type']})\n")
                    f.write(f"Predicted preferred: {row['pred_preferred']} ({row['pred_type']})\n")
                    f.write(f"Energy difference: {row['energy_diff']:.2f} eV\n")
                    f.write(f"DFT all preferred: {row['dft_all_preferred']}\n")
                    f.write(f"Predicted all preferred: {row['pred_all_preferred']}\n")
        
        # Create coordination comparison plots
        plot_coordination_comparison(coord_results, root, output_suffix, args.model, args.Y, row_str)
        print("Coordination preference analysis completed")
    else:
        print(f"{YELLOW}No valid results for coordination preference analysis{ENDC}")

    print_time("Results saving completed", time.time() - save_start)
    print_time("Total execution time", time.time() - start_time)

if __name__ == '__main__':
    main() 