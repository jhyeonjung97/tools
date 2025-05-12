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
from sklearn.impute import SimpleImputer
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
import xgboost as xgb
import lightgbm as lgb

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
    'form': 'Formation Energy (eV)',
    'coh': 'Cohesive Energy (eV)',
}

def print_time(message, time_value):
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

# Drop-column importance 구현
from sklearn.base import clone
from sklearn.metrics import mean_absolute_error

def feature_importance_drop_column(X_train, y_train, X_test, y_test, model, feature_names, model_type='gpr'):
    """
    Drop-column importance: 각 feature를 제외하고 모델을 다시 학습하여 test MAE가 얼마나 증가하는지로 중요도를 계산
    """
    # 전체 feature로 학습한 모델의 test MAE
    model_full = clone(model)
    model_full.fit(X_train, y_train)
    base_mae = mean_absolute_error(y_test, model_full.predict(X_test))
    importances = []
    for i, feat in enumerate(feature_names):
        drop_feats = [f for j, f in enumerate(feature_names) if j != i]
        X_train_drop = X_train[drop_feats]
        X_test_drop = X_test[drop_feats]
        model_drop = clone(model)
        model_drop.fit(X_train_drop, y_train)
        drop_mae = mean_absolute_error(y_test, model_drop.predict(X_test_drop))
        importances.append(drop_mae - base_mae)
    return pd.Series(importances, index=feature_names)

# 이하 bulk_pred_cfse.py의 main 함수 및 기타 함수는 그대로 복사, feature_importance 부분만 위 함수로 대체

def main():
    start_time = time.time()
    print("Starting bulk prediction analysis (drop-column importance)...")
    
    parser = argparse.ArgumentParser(description='Bulk prediction using GPR, GBR, RF, LR, XGBoost, or LightGBM (drop-column importance)')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr', 'xgb', 'lgb'], default='gpr',
                      help='Model type to use (gpr, gbr, rf, lr, xgb, or lgb)')
    parser.add_argument('--Y', type=str, choices=['form', 'coh'], default='form',
                      help='Target column from bulk_data_total.csv (form or coh)')
    parser.add_argument('--X', nargs='+', default=[
        'OS', 'CN', 'numb', 'chg', 'mag', 'volume', 'l_bond', 'madelung', 
        'ICOHP', 'ICOHPo', 'ICOHPc', 'ICOBI', 'ICOBIo', 'ICOBIc', 'ICOOP', 'ICOOPo', 'ICOOPc', 
        'ion-1', 'ion', 'ion+1', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1', 
        'pauling', 'Natom', 'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform', 'n_electrons', 'd_electrons',
        'base_cfse', 'ee_repulsion', 'jt_effect', 'field_strength', 'cfse', 'exchange_stabilization'
    ], help='List of feature columns from bulk_data_total.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size (default: 0.2)')
    parser.add_argument('--random_state', type=int, default=50, help='Random state for reproducibility (default: 50)')
    parser.add_argument('--energy_threshold', type=float, default=0.2,
                       help='Energy threshold (eV) for considering multiple coordinations in preference analysis')
    parser.add_argument('--corr_threshold', type=float, default=1.0,
                       help='Correlation threshold for feature selection (0.7, 0.8, 0.9 등)')
    parser.add_argument('--save', action='store_true', help='Save results to file')
    args = parser.parse_args()
    
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    
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
    
    output_suffix = args.output
    row_str = ''.join(sorted(args.row)) if args.row else 'all'
    threshold_str = str(int(args.corr_threshold * 100)) if args.corr_threshold != 1.0 else '00'
    log_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}.log')
    tsv_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}.tsv')
    png_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}.png')
    importance_png_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}_importance.png')
    json_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}.json')

    try:
        print("Loading data...")
        load_start = time.time()
        df_bulk = pd.read_csv(os.path.join(root, 'bulk_data_cfse.csv'), index_col=0)
        print_time("Data loading completed", time.time() - load_start)
    except FileNotFoundError as e:
        print(f"{RED}Error: Required data file not found: {e}{ENDC}")
        exit(1)

    print("Preprocessing data...")
    preprocess_start = time.time()
    df = df_bulk.copy()

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]

    X = df[args.X].astype(float)
    y = df[args.Y].astype(float)
    valid_indices = ~y.isna()
    X = X[valid_indices]
    y = y[valid_indices]
    df = df[valid_indices]

    imputer = SimpleImputer(strategy='median')
    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns, index=X.index)

    print("Calculating correlation matrix...")
    corr_matrix = df[args.X + [args.Y]].corr()
    print("Selecting features based on correlation...")
    selected_features = args.X  # drop-column importance는 feature selection과 무관하게 전체 feature 사용

    X = df[selected_features].astype(float)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.random_state
    )

    print(f"\nTraining and test dataset sizes:")
    print(f"Training data: {len(X_train)} samples")
    print(f"Test data: {len(X_test)} samples")
    print(f"Total data: {len(X_train) + len(X_test)} samples")

    print("Scaling features...")
    scale_start = time.time()
    scaler = StandardScaler()
    X_train_scaled = pd.DataFrame(
        scaler.fit_transform(X_train),
        columns=X_train.columns,
        index=X_train.index
    )
    X_test_scaled = pd.DataFrame(
        scaler.transform(X_test),
        columns=X_test.columns,
        index=X_test.index
    )
    print_time("Feature scaling completed", time.time() - scale_start)

    print(f"Setting up {args.model.upper()} model...")
    model_start = time.time()
    if args.model == 'gpr':
        kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e5)) * RationalQuadratic(
            length_scale=1.0,
            alpha=1.0,
            length_scale_bounds=(1e-3, 1e6),
            alpha_bounds=(1e-3, 1e6)
        ) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-3, 1e5))
        model = MyGPR(
            kernel=kernel,
            random_state=args.random_state,
            n_restarts_optimizer=20,
            alpha=1e-3,
            optimizer='fmin_l_bfgs_b',
            normalize_y=True,
            max_iter=1000,
            gtol=1e-3
        )
    elif args.model == 'gbr':
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
    elif args.model == 'rf':
        model = RandomForestRegressor(
            n_estimators=100,
            random_state=args.random_state
        )
    elif args.model == 'lr':
        model = LinearRegression()
    elif args.model == 'xgb':
        model = xgb.XGBRegressor(
            objective='reg:squarederror',
            random_state=args.random_state,
            n_estimators=200,
            learning_rate=0.05,
            max_depth=6,
            min_child_weight=1,
            subsample=0.8,
            colsample_bytree=0.8,
            gamma=0,
            reg_alpha=0,
            reg_lambda=1
        )
    else:  # lgb
        model = lgb.LGBMRegressor(
            n_estimators=1000,
            learning_rate=0.01,
            num_leaves=31,
            random_state=args.random_state,
            n_jobs=1,
            verbose=-1
        )

    print("Performing cross-validation...")
    cv_start = time.time()
    kf = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    X_train_scaled_values = X_train_scaled.values
    cv_r2 = cross_val_score(model, X_train_scaled_values, y_train, cv=kf, scoring='r2')
    cv_mae = -cross_val_score(model, X_train_scaled_values, y_train, cv=kf, scoring='neg_mean_absolute_error')
    cv_mse = -cross_val_score(model, X_train_scaled_values, y_train, cv=kf, scoring='neg_mean_squared_error')
    print_time("Cross-validation completed", time.time() - cv_start)
    print(f"{MAGENTA}Cross-validation R2 scores: {cv_r2}{ENDC}")
    print(f"{MAGENTA}Mean CV R2 score: {cv_r2.mean():.4f} (+/- {cv_r2.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation MAE scores: {cv_mae}{ENDC}")
    print(f"{MAGENTA}Mean CV MAE score: {cv_mae.mean():.4f} (+/- {cv_mae.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation MSE scores: {cv_mse}{ENDC}")
    print(f"{MAGENTA}Mean CV MSE score: {cv_mse.mean():.4f} (+/- {cv_mse.std() * 2:.4f}){ENDC}")

    print("Fitting model...")
    model.fit(X_train_scaled, y_train)
    print("Making predictions...")
    pred_start = time.time()
    y_pred_train = model.predict(X_train_scaled)
    y_pred_test = model.predict(X_test_scaled)
    print_time("Predictions completed", time.time() - pred_start)

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

    print("Calculating drop-column feature importance...")
    importance_start = time.time()
    importance = feature_importance_drop_column(X_train_scaled, y_train, X_test_scaled, y_test, model, selected_features, args.model)
    print_time("Feature importance calculation completed", time.time() - importance_start)

    print("Saving results...")
    save_start = time.time()
    with open(log_path, 'w') as f:
        f.write("Training Metrics:\n")
        f.write(f"train R2: {metrics['train']['r2']:.4f}\n")
        f.write(f"train MAE: {metrics['train']['mae']:.4f}\n")
        f.write(f"train MSE: {metrics['train']['mse']:.4f}\n")
        f.write(f"train RMSE: {metrics['train']['rmse']:.4f}\n")
        f.write(f"train Max Error: {metrics['train']['max_error']:.4f}\n\n")
        f.write("Test Metrics:\n")
        f.write(f"test R2: {metrics['test']['r2']:.4f}\n")
        f.write(f"test MAE: {metrics['test']['mae']:.4f}\n")
        f.write(f"test MSE: {metrics['test']['mse']:.4f}\n")
        f.write(f"test RMSE: {metrics['test']['rmse']:.4f}\n")
        f.write(f"test Max Error: {metrics['test']['max_error']:.4f}\n\n")
        f.write("Cross-validation Results:\n")
        f.write(f"R2 scores: {cv_r2}\n")
        f.write(f"Mean R2 score: {cv_r2.mean():.4f} (+/- {cv_r2.std() * 2:.4f})\n")
        f.write(f"MAE scores: {cv_mae}\n")
        f.write(f"Mean MAE score: {cv_mae.mean():.4f} (+/- {cv_mae.std() * 2:.4f})\n")
        f.write(f"MSE scores: {cv_mse}\n")
        f.write(f"Mean MSE score: {cv_mse.mean():.4f} (+/- {cv_mse.std() * 2:.4f})\n")
        f.write("\nDrop-column Feature Importance (test MAE increase):\n")
        for feat, imp in importance.items():
            f.write(f"{feat}: {imp:.4f}\n")
    print(f"{BLUE}Log file saved as {log_path}{ENDC}")

    output_path = os.path.join(root, f'{args.Y}_pred_cfse_{args.model}{threshold_str}_{row_str}_{output_suffix}_importance.png')
    plt.figure(figsize=(10, 4))
    importance.sort_values(ascending=False).plot(kind='bar')
    plt.ylabel('Test MAE Increase (Drop-Column)')
    plt.title('Drop-Column Feature Importance')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"{BLUE}Feature importance plot saved as {output_path}{ENDC}")

if __name__ == '__main__':
    print("Script started...")
    print("Importing required libraries...")
    main() 