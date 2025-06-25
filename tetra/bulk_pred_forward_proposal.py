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
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, max_error
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.utils.optimize import _check_optimize_result
import scipy.optimize
import socket
import xgboost as xgb
import lightgbm as lgb
from sklearn.base import clone
import matplotlib.ticker as ticker

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

def forward_feature_selection(X_train, y_train, X_test, y_test, model, feature_names):
    """
    Sequential Forward Selection: feature를 하나씩 추가하며 test MAE가 가장 낮아지는 feature를 선택
    """
    selected = []
    remaining = list(feature_names)
    history = []
    for step in range(min(5, len(feature_names))):  # 최대 5개 feature까지만 선택
        best_test_mae = np.inf
        best_test_mse = np.inf
        best_feat = None
        for feat in remaining:
            candidate = selected + [feat]
            model_ = clone(model)
            model_.fit(X_train[candidate], y_train)
            y_pred_test = model_.predict(X_test[candidate])
            test_mae = mean_absolute_error(y_test, y_pred_test)
            test_mse = mean_squared_error(y_test, y_pred_test)
            if test_mae < best_test_mae:
                best_test_mae = test_mae
                best_test_mse = test_mse
                best_feat = feat
        selected.append(best_feat)
        remaining.remove(best_feat)
        history.append({
            'n_features': len(selected),
            'features': selected.copy(),
            'test_mae': best_test_mae,
            'test_mse': best_test_mse
        })
        print(f"Step {len(selected)}: Added '{best_feat}', Test MAE: {best_test_mae:.4f}, Test MSE: {best_test_mse:.4f}")
    return history

def main():
    start_time = time.time()
    print("Starting bulk prediction analysis (forward selection)...")
    parser = argparse.ArgumentParser(description='Bulk prediction with forward feature selection')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr', 'xgb', 'lgb'], default='gpr', help='Model type to use')
    parser.add_argument('--Y', type=str, choices=['form', 'coh'], default='form', help='Target column')
    parser.add_argument('--X', nargs='+', required=True, help='List of feature columns')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination')
    parser.add_argument('--output', type=str, default='forward', help='Output filename suffix')
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
    log_path = os.path.join(root, f'{args.Y}_forward_{args.model}_{row_str}_{output_suffix}.log')
    png_path = os.path.join(root, f'{args.Y}_forward_{args.model}_{row_str}_{output_suffix}.png')

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
    print("Running forward feature selection...")
    history = forward_feature_selection(X_train_scaled, y_train, X_test_scaled, y_test, model, args.X)
    # 로그 저장
    with open(log_path, 'w') as f:
        for record in history:
            f.write(f"n_features: {record['n_features']}, features: {record['features']}, test_mae: {record['test_mae']:.4f}, test_mse: {record['test_mse']:.4f}\n")
    print(f"{BLUE}Forward selection log saved as {log_path}{ENDC}")
    # MAE/MSE plot
    test_maes = [h['test_mae'] for h in history]
    test_rmses = [np.sqrt(h['test_mse']) for h in history]  # MSE의 제곱근을 취해 RMSE 계산
    feature_names_order = [h['features'][-1] for h in history]
    fig, ax = plt.subplots(figsize=(4, 2.7))
    ax.plot(range(1, 6), test_maes[:5], marker='o', label='Test MAE')
    ax.plot(range(1, 6), test_rmses[:5], marker='s', label='Test RMSE')
    ax.set_xlabel('Feature Added')
    ax.set_ylabel('Error (eV)')
    ax.set_ylim(0.2, 1.5)
    # ax.set_xticks(range(1, 6))
    ax.set_xticklabels([])  # 아래쪽 x축 숫자 숨김
    # ax.tick_params(axis='x', bottom=False)  # 아래쪽 x축 눈금선만 숨김
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    # 위쪽에 feature 이름 표시
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(1, 6))
    feature_names_order = ['outer\nelectron', 'ionization\nenergy', 'evaporation\nenthalpy', 'normalized\nICOHP', 'Bader\ncharge']
    ax2.set_xticklabels(feature_names_order, rotation=45)
    ax.legend()
    fig.tight_layout()
    plt.savefig(png_path, transparent=True, dpi=300)
    plt.close()
    print(f"{BLUE}Forward selection MAE/MSE plot saved as {png_path}{ENDC}")

if __name__ == '__main__':
    main() 
