import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic, WhiteKernel, ConstantKernel
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.utils.optimize import _check_optimize_result
import scipy.optimize
import socket
import xgboost as xgb
import lightgbm as lgb
from sklearn.base import clone

RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'
CYAN = '\033[96m'
ENDC = '\033[0m'
BOLD = '\033[1m'

os.environ['JOBLIB_START_METHOD'] = 'fork'

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

def drop_column_importance(X_train, y_train, X_test, y_test, model, feature_names):
    """
    Drop-column importance: 각 feature를 하나씩 제외하고 모델을 학습하여 test MAE, test MSE가 얼마나 증가하는지로 중요도를 계산
    """
    # 전체 feature로 학습한 모델의 test MAE, test MSE
    model_full = clone(model)
    model_full.fit(X_train, y_train)
    base_mae = mean_absolute_error(y_test, model_full.predict(X_test))
    base_mse = mean_squared_error(y_test, model_full.predict(X_test))
    importances = []
    for i, feat in enumerate(feature_names):
        drop_feats = [f for j, f in enumerate(feature_names) if j != i]
        X_train_drop = X_train[drop_feats]
        X_test_drop = X_test[drop_feats]
        model_drop = clone(model)
        model_drop.fit(X_train_drop, y_train)
        drop_mae = mean_absolute_error(y_test, model_drop.predict(X_test_drop))
        drop_mse = mean_squared_error(y_test, model_drop.predict(X_test_drop))
        importances.append({
            'feature': feat,
            'mae_increase': drop_mae - base_mae,
            'mse_increase': drop_mse - base_mse,
            'drop_mae': drop_mae,
            'drop_mse': drop_mse
        })
        print(f"Drop '{feat}': Test MAE increase: {drop_mae - base_mae:.4f}, Test MSE increase: {drop_mse - base_mse:.4f}")
    return importances, base_mae, base_mse

def main():
    start_time = time.time()
    print("Starting bulk prediction analysis (drop-column importance)...")
    parser = argparse.ArgumentParser(description='Bulk prediction with drop-column feature importance')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr', 'xgb', 'lgb'], default='gpr', help='Model type to use')
    parser.add_argument('--Y', type=str, choices=['form', 'coh'], default='form', help='Target column')
    parser.add_argument('--X', nargs='+', required=True, help='List of feature columns')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination')
    parser.add_argument('--output', type=str, default='drop', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size')
    parser.add_argument('--random_state', type=int, default=50, help='Random state')
    args = parser.parse_args()

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
    log_path = os.path.join(root, f'{args.Y}_drop_{args.model}_{row_str}_{output_suffix}.log')
    png_path = os.path.join(root, f'{args.Y}_drop_{args.model}_{row_str}_{output_suffix}.png')

    print("Loading data...")
    df_bulk = pd.read_csv(os.path.join(root, 'bulk_data_cfse.csv'), index_col=0)
    # feature 이름 전처리 추가
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
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
    X_train, X_test, y_train, y_test = train_test_split(
        X_imputed, y, test_size=args.test_size, random_state=args.random_state
    )
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
    print(f"Training data: {len(X_train_scaled)} samples, Test data: {len(X_test_scaled)} samples")
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
            n_restarts_optimizer=10,
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
    print("Running drop-column feature importance...")
    importances, base_mae, base_mse = drop_column_importance(X_train_scaled, y_train, X_test_scaled, y_test, model, args.X)
    # 로그 저장
    with open(log_path, 'w') as f:
        f.write(f"Base test MAE: {base_mae:.4f}, Base test MSE: {base_mse:.4f}\n")
        for record in importances:
            f.write(f"feature: {record['feature']}, mae_increase: {record['mae_increase']:.4f}, mse_increase: {record['mse_increase']:.4f}, drop_mae: {record['drop_mae']:.4f}, drop_mse: {record['drop_mse']:.4f}\n")
    print(f"{BLUE}Drop-column importance log saved as {log_path}{ENDC}")
    # MAE/MSE plot (막대바, 겹치지 않게)
    features = [r['feature'] for r in importances]
    mae_increases = [r['mae_increase'] for r in importances]
    mse_increases = [r['mse_increase'] for r in importances]
    # 중요도(오름차순)로 정렬 (값이 클수록 중요)
    sorted_idx = np.argsort(mae_increases)[::-1]
    features_sorted = [features[i] for i in sorted_idx]
    mae_increases_sorted = [mae_increases[i] for i in sorted_idx]
    mse_increases_sorted = [mse_increases[i] for i in sorted_idx]
    x = np.arange(len(features_sorted))
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 4))
    rects1 = ax.bar(x - width/2, mae_increases_sorted, width, label='Test MAE Increase')
    rects2 = ax.bar(x + width/2, mse_increases_sorted, width, label='Test MSE Increase')
    ax.set_xlabel('Feature')
    ax.set_ylabel('Error Increase')
    ax.set_xticks(x)
    ax.set_xticklabels(features_sorted, rotation=90)
    ax.legend()
    fig.tight_layout()
    plt.savefig(png_path)
    plt.close()
    print(f"{BLUE}Drop-column importance plot saved as {png_path}{ENDC}")

if __name__ == '__main__':
    main() 