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
from sklearn.gaussian_process.kernels import RationalQuadratic, WhiteKernel, ConstantKernel
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, max_error
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.utils.optimize import _check_optimize_result
import scipy.optimize
import socket

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
    'grosspop': 'Gross Population (Loewdin e⁻)',
}

def print_time(message, time_value):
    """Print time-related information with color coding"""
    print(f"{BOLD}{message}: {time_value:.2f} seconds{ENDC}")

class MyGPR(GaussianProcessRegressor):
    def __init__(self, kernel=None, alpha=1e-10, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=0, normalize_y=False, copy_X_train=True, random_state=None, max_iter=1e06, gtol=1e-05):
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

def feature_importance(X, y, model, feature_names, model_type='gpr'):
    """Calculate feature importance"""
    if model_type == 'gpr':
        importance = []
        for i in tqdm(range(X.shape[1]), desc="Calculating feature importance"):
            X_permuted = X.copy()
            X_permuted[:, i] = np.random.permutation(X_permuted[:, i])
            score = model.score(X_permuted, y)
            importance.append(1 - score)
        return pd.Series(importance, index=feature_names)
    elif model_type in ['gbr', 'rf']:
        return pd.Series(model.feature_importances_, index=feature_names)
    else:  # lr
        return pd.Series(np.abs(model.coef_), index=feature_names)

def get_prediction_std(model, X, n_iterations=100, model_type='gpr'):
    """
    Estimate prediction uncertainty
    Args:
        model: Trained model (GPR, GBR, RF, or LR)
        X: Input features
        n_iterations: Number of bootstrap iterations (for GBR/RF)
        model_type: Type of model ('gpr', 'gbr', 'rf', or 'lr')
    Returns:
        Standard deviation of predictions
    """
    if model_type == 'gpr':
        _, std = model.predict(X, return_std=True)
        return std
    elif model_type in ['gbr', 'rf']:
        predictions = []
        n_estimators = len(model.estimators_)
        
        for i in range(n_iterations):
            indices = np.random.choice(n_estimators, size=n_estimators//2, replace=False)
            pred = np.zeros(X.shape[0])
            for idx in indices:
                pred += model.estimators_[idx].predict(X)
            pred /= len(indices)
            predictions.append(pred)
        return np.std(predictions, axis=0)
    else:  # lr
        # For LR, we'll use the standard deviation of the residuals as a simple uncertainty estimate
        y_pred = model.predict(X)
        residuals = y_pred - model.predict(X)
        return np.std(residuals) * np.ones(X.shape[0])

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
            model_selected = MyGPR(kernel=kernel, random_state=42)
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
            model_selected = MyGPR(kernel=kernel, random_state=42)
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

def plot_training_metrics(train_mae, train_mse, test_mae, test_mse, output_dir, output_suffix, model_type):
    """Plot training and test MAE and MSE over iterations in separate files."""
    # Get Y and row information from output_suffix
    parts = output_suffix.split('_')
    Y = parts[0]  # Y는 output_suffix의 첫 번째 부분
    row_str = parts[2]  # row 정보는 세 번째 부분
    
    # Plot MAE
    plt.figure(figsize=(5, 4))
    plt.plot(train_mae, label='Train MAE', color='blue')
    plt.plot(test_mae, label='Test MAE', color='orange')
    plt.xlabel('Iterations')
    plt.ylabel('MAE')
    plt.legend()
    plt.tight_layout()
    mae_path = os.path.join(output_dir, f'{Y}_pred_{model_type}_{row_str}_mae.png')
    plt.savefig(mae_path)
    plt.close()
    print(f"{BLUE}MAE training metrics plot saved as {mae_path}{ENDC}")
    
    # Plot MSE
    plt.figure(figsize=(5, 4))
    plt.plot(train_mse, label='Train MSE', color='blue')
    plt.plot(test_mse, label='Test MSE', color='orange')
    plt.xlabel('Iterations')
    plt.ylabel('MSE')
    plt.legend()
    plt.tight_layout()
    mse_path = os.path.join(output_dir, f'{Y}_pred_{model_type}_{row_str}_mse.png')
    plt.savefig(mse_path)
    plt.close()
    print(f"{BLUE}MSE training metrics plot saved as {mse_path}{ENDC}")

def main():
    start_time = time.time()
    print("Starting bulk prediction analysis...")
    
    parser = argparse.ArgumentParser(description='Bulk prediction using GPR, GBR, RF, or LR')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr'], default='gpr',
                      help='Model type to use (gpr, gbr, rf, or lr)')
    parser.add_argument('--Y', type=str, choices=['form', 'coh'], default='form',
                      help='Target column from bulk_data.csv (form or coh)')
    parser.add_argument('--X', nargs='+', default=[
        'numb', 'chg', 'mag', 'volume', 'l_bond', 'n_bond',
        'grosspop', 'madelung', 'ICOHPm', 'ICOHPn', 'ICOBIm', 'ICOBIn', 'ICOOPm', 'ICOOPn', 
        'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
        'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size (default: 0.2)')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility (default: 42)')
    parser.add_argument('--threshold', type=float, default=0.2,
                       help='Energy threshold (eV) for considering multiple coordinations in preference analysis')
    args = parser.parse_args()
    
    # Convert feature names if they start with ICOHP or ICOOP (prepend '-')
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    
    # 서버 주소 가져오기
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
    
    # Define paths for all output files
    output_suffix = args.output
    row_str = ''.join(sorted(args.row)) if args.row else 'all'  # row 정보를 문자열로 변환
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
    df = pd.merge(df_bulk, df_mend, left_on='metal', right_index=True, suffixes=('_bulk', '_mend'))
    df = df.rename(columns={'row_bulk': 'row', 'numb_mend': 'numb'})
    df = df.drop(columns=['row_mend', 'numb_bulk'])
    df = df[df['row'] != 'fm']

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    # Drop rows with NaN in any relevant column
    all_columns = args.X + [args.Y]
    df = df.dropna(subset=all_columns)

    X = df[args.X].astype(float)
    y = df[args.Y].astype(float)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=args.test_size, random_state=args.random_state)

    # Scale the features
    print("Scaling features...")
    scale_start = time.time()
    scaler = RobustScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    print_time("Feature scaling completed", time.time() - scale_start)

    # Define and fit model
    print(f"Setting up {args.model.upper()} model...")
    model_start = time.time()
    
    if args.model == 'gpr':
        # GPR model setup
        kernel = ConstantKernel(1.0, constant_value_bounds=(1e-5, 1e5)) * RationalQuadratic(
            length_scale=1.0,
            alpha=1.0,
            length_scale_bounds=(1e-5, 1e5),
            alpha_bounds=(1e-5, 1e5)
        ) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-5, 1e5))
        
        model = MyGPR(
            kernel=kernel,
            random_state=args.random_state,
            n_restarts_optimizer=50,
            alpha=1e-5,
            optimizer='fmin_l_bfgs_b',
            normalize_y=True,
            max_iter=1e07,
            gtol=1e-5
        )

        # Bayesian Optimization for GPR
        search_space = {
            'kernel__k1__k2__length_scale': Real(1e-5, 1e5, prior='log-uniform'),
            'kernel__k1__k2__alpha': Real(1e-5, 1e5, prior='log-uniform'),
            'kernel__k2__noise_level': Real(1e-5, 1e5, prior='log-uniform')
        }
    elif args.model == 'gbr':
        # GBR model setup
        model = GradientBoostingRegressor(
            n_estimators=200,
            learning_rate=0.05,
            max_depth=4,
            min_samples_split=10,
            min_samples_leaf=5,
            subsample=0.85,
            max_features=0.8,
            validation_fraction=0.15,
            n_iter_no_change=15,
            tol=1e-5,
            random_state=args.random_state
        )

        # Bayesian Optimization for GBR
        search_space = {
            'n_estimators': Integer(150, 300),
            'learning_rate': Real(0.03, 0.08),
            'max_depth': Integer(3, 5),
            'min_samples_split': Integer(8, 15),
            'min_samples_leaf': Integer(4, 8),
            'subsample': Real(0.8, 0.9),
            'max_features': Real(0.7, 0.9),
        }
    elif args.model == 'rf':
        # RF model setup
        model = RandomForestRegressor(
            n_estimators=100,
            random_state=args.random_state
        )

        # Bayesian Optimization for RF
        search_space = {
            'n_estimators': Integer(50, 500),
            'max_depth': Integer(1, 20),
            'min_samples_split': Integer(2, 20),
            'min_samples_leaf': Integer(1, 10),
        }
    else:  # lr
        # LR model setup
        model = LinearRegression()
        # No hyperparameter optimization needed for LR
        search_space = {}

    # Perform Bayesian Optimization
    print("Performing Bayesian Optimization...")
    bayes_start = time.time()
    
    if args.model != 'lr':  # LR 모델이 아닐 때만 Bayesian Optimization 수행
        cv = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
        bayes_search = BayesSearchCV(
            model,
            search_space,
            n_iter=50,
            cv=cv,
            scoring='neg_mean_squared_error',
            n_jobs=1,
            random_state=args.random_state,
            verbose=1
        )
        bayes_search.fit(X_train_scaled, y_train)
        print_time("Bayesian Optimization completed", time.time() - bayes_start)
        print(f"{MAGENTA}Best parameters: {bayes_search.best_params_}{ENDC}")
        
        # Use best model
        model = bayes_search.best_estimator_
    else:  # LR 모델인 경우
        model.fit(X_train_scaled, y_train)
        bayes_search = None
        print_time("Model fitting completed", time.time() - bayes_start)

    # Perform cross-validation
    print("Performing cross-validation...")
    cv_start = time.time()
    kf = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    
    # Calculate multiple metrics for cross-validation
    cv_r2 = cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='r2')
    cv_mae = -cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='neg_mean_absolute_error')
    cv_mse = -cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='neg_mean_squared_error')
    
    print_time("Cross-validation completed", time.time() - cv_start)
    print(f"{MAGENTA}Cross-validation R2 scores: {cv_r2}{ENDC}")
    print(f"{MAGENTA}Mean CV R2 score: {cv_r2.mean():.4f} (+/- {cv_r2.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation MAE scores: {cv_mae}{ENDC}")
    print(f"{MAGENTA}Mean CV MAE score: {cv_mae.mean():.4f} (+/- {cv_mae.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation MSE scores: {cv_mse}{ENDC}")
    print(f"{MAGENTA}Mean CV MSE score: {cv_mse.mean():.4f} (+/- {cv_mse.std() * 2:.4f}){ENDC}")
    
    # Make predictions
    print("Making predictions...")
    pred_start = time.time()
    if args.model == 'gpr':
        y_pred_train, std_train = model.predict(X_train_scaled, return_std=True)
        y_pred_test, std_test = model.predict(X_test_scaled, return_std=True)
    elif args.model == 'gbr':
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        std_train = get_prediction_std(model, X_train_scaled, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled, model_type=args.model)
    elif args.model == 'rf':
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        std_train = get_prediction_std(model, X_train_scaled, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled, model_type=args.model)
    else:  # lr
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        std_train = get_prediction_std(model, X_train_scaled, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled, model_type=args.model)
    print_time("Predictions completed", time.time() - pred_start)

    # Calculate metrics
    print("Calculating metrics...")
    metrics_start = time.time()
    metrics = {
        'train': {
            'r2': r2_score(y_train, y_pred_train),
            'mae': mean_absolute_error(y_train, y_pred_train),
            'mse': mean_squared_error(y_train, y_pred_train),
            'max_error': max_error(y_train, y_pred_train),
            'rmse': np.sqrt(mean_squared_error(y_train, y_pred_train))
        },
        'test': {
            'r2': r2_score(y_test, y_pred_test),
            'mae': mean_absolute_error(y_test, y_pred_test),
            'mse': mean_squared_error(y_test, y_pred_test),
            'max_error': max_error(y_test, y_pred_test),
            'rmse': np.sqrt(mean_squared_error(y_test, y_pred_test))
        },
        'cv': {
            'mean': cv_r2.mean(),
            'std': cv_r2.std(),
            'scores': cv_r2.tolist()
        }
    }
    print_time("Metrics calculation completed", time.time() - metrics_start)

    # Calculate feature importance
    print("Calculating feature importance...")
    importance_start = time.time()
    importance = feature_importance(X_train_scaled, y_train, model, args.X, args.model)
    cumulative_scores = plot_feature_importance(X_train_scaled, X_test_scaled, y_train, y_test, args.X, importance_png_path, args.model)
    print_time("Feature importance calculation completed", time.time() - importance_start)

    # Save results
    print("Saving results...")
    save_start = time.time()

    # Save metrics and cross-validation results
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
        f.write(f"Max Error: {metrics['test']['max_error']:.4f}\n\n")
        f.write("Cross-validation Results:\n")
        f.write(f"R2 scores: {cv_r2}\n")
        f.write(f"Mean R2 score: {cv_r2.mean():.4f} (+/- {cv_r2.std() * 2:.4f})\n")
        f.write(f"MAE scores: {cv_mae}\n")
        f.write(f"Mean MAE score: {cv_mae.mean():.4f} (+/- {cv_mae.std() * 2:.4f})\n")
        f.write(f"MSE scores: {cv_mse}\n")
        f.write(f"Mean MSE score: {cv_mse.mean():.4f} (+/- {cv_mse.std() * 2:.4f})\n")
        
        # Add Feature Importance Analysis
        f.write("\nFeature Importance Analysis:\n")
        for score in cumulative_scores:
            f.write(f"\nNumber of features: {score['n_features']}\n")
            f.write(f"Features used: {', '.join(score['features'])}\n")
            f.write(f"Train MAE: {score['mae_train']:.4f}\n")
            f.write(f"Test MAE: {score['mae_test']:.4f}\n")
    print(f"{BLUE}Log file saved as {log_path}{ENDC}")

    # Save predictions
    df_result = pd.DataFrame({
        'metal': df.index,
        'row': df['row'],
        'coord': df['coord'],
        'Y_true': y,
        'Y_pred': model.predict(scaler.transform(X)),
        'std': get_prediction_std(model, scaler.transform(X), model_type=args.model)
    })
    
    # Extract actual metal name from the index
    df_result['metal'] = df_result.apply(lambda row: df.loc[row['metal'], 'metal'], axis=1)
    
    df_result.to_csv(tsv_path, sep='\t', index=False)
    print(f"{BLUE}TSV file saved as {tsv_path}{ENDC}")

    # Plot parity plot with uncertainty
    print("Creating parity plot...")
    plt.figure(figsize=(10, 8))
    
    # Create DataFrames for train and test sets
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
    
    # Define color and marker mappings
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}
    coord_map = {'WZ': '>', 'ZB': '<', 'TN': 'o', 'PD': 'o', 'NB': 's', 'RS': 'd', 'LT': 'h'}
    
    # Plot training data with error bars
    for r in df['row'].unique():
        for c in df['coord'].unique():
            # Plot training data
            subset_train = df_train[(df_train['row'] == r) & (df_train['coord'] == c)]
            if not subset_train.empty:
                row_features = subset_train[args.X].astype(float)
                row_features_scaled = scaler.transform(row_features)
                if args.model == 'gpr':
                    y_pred_train, std_train = model.predict(row_features_scaled, return_std=True)
                else:
                    y_pred_train = model.predict(row_features_scaled)
                    std_train = get_prediction_std(model, row_features_scaled, model_type=args.model)
                
                plt.errorbar(
                    subset_train['Y'].values,
                    y_pred_train,
                    yerr=std_train,
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
                               (row_data['Y'], y_pred_train[i]), 
                               fontsize=8)
            
            # Plot test data
            subset_test = df_test[(df_test['row'] == r) & (df_test['coord'] == c)]
            if not subset_test.empty:
                row_features = subset_test[args.X].astype(float)
                row_features_scaled = scaler.transform(row_features)
                if args.model == 'gpr':
                    y_pred_test, std_test = model.predict(row_features_scaled, return_std=True)
                else:
                    y_pred_test = model.predict(row_features_scaled)
                    std_test = get_prediction_std(model, row_features_scaled, model_type=args.model)
                
                plt.errorbar(
                    subset_test['Y'].values,
                    y_pred_test,
                    yerr=std_test,
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
                               (row_data['Y'], y_pred_test[i]), 
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

    # Save all results in JSON format
    results = {
        'metrics': metrics,
        'model_params': {
            'type': args.model,
            'best_params': bayes_search.best_params_ if bayes_search else None
        },
        'feature_names': args.X,
        'target': args.Y,
        'feature_importance': importance.to_dict()
    }
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"{BLUE}Results saved as {json_path}{ENDC}")

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

    # Calculate training metrics over iterations for GBR
    if args.model == 'gbr':
        print("Calculating training metrics over iterations...")
        train_mae_list = []
        train_mse_list = []
        test_mae_list = []
        test_mse_list = []

        for y_pred_train in model.staged_predict(X_train_scaled):
            # 학습 데이터에 대한 메트릭 계산
            train_mae = mean_absolute_error(y_train, y_pred_train)
            train_mse = mean_squared_error(y_train, y_pred_train)
            
            train_mae_list.append(train_mae)
            train_mse_list.append(train_mse)
            
            # 테스트 데이터에 대한 메트릭 계산
            y_pred_test = model.predict(X_test_scaled)
            test_mae = mean_absolute_error(y_test, y_pred_test)
            test_mse = mean_squared_error(y_test, y_pred_test)
            
            test_mae_list.append(test_mae)
            test_mse_list.append(test_mse)

        # Plot training metrics
        plot_training_metrics(train_mae_list, train_mse_list, test_mae_list, test_mse_list, root, f'{args.Y}_pred_{args.model}_{row_str}_{output_suffix}', args.model)

    # Print total execution time
    total_time = time.time() - start_time
    print("Execution Summary:")
    print_time("Total execution time", total_time)
    print_time("Model training and evaluation time", time.time() - model_start)
    
    # Save execution time to log file
    with open(log_path, 'a') as f:
        f.write("\nExecution Times:\n")
        f.write(f"Total execution time: {total_time:.2f} seconds\n")
        f.write(f"Data loading: {time.time() - load_start:.2f} seconds\n")
        f.write(f"Data preprocessing: {time.time() - preprocess_start:.2f} seconds\n")
        f.write(f"Feature scaling: {time.time() - scale_start:.2f} seconds\n")
        if args.model in ['gpr', 'gbr', 'rf']:
            f.write(f"Bayesian Optimization: {time.time() - bayes_start:.2f} seconds\n")
        f.write(f"Cross-validation: {time.time() - cv_start:.2f} seconds\n")
        f.write(f"Predictions: {time.time() - pred_start:.2f} seconds\n")
        f.write(f"Metrics calculation: {time.time() - metrics_start:.2f} seconds\n")
        f.write(f"Feature importance calculation: {time.time() - importance_start:.2f} seconds\n")
        f.write(f"Results saving and plotting: {time.time() - save_start:.2f} seconds\n")

if __name__ == '__main__':
    main() 