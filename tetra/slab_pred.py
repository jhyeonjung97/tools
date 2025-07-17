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

def main():
    start_time = time.time()
    print("Starting surface adsorption prediction analysis...")
    hostname = socket.gethostname()
    user_name = os.getlogin()
    if hostname == 'PC102616':
        root = '/Users/jiuy97/Desktop/8_V_slab/figures'
    elif user_name == 'jiuy97':
        root = '/pscratch/sd/j/jiuy97/8_V_slab/figures'
    elif user_name == 'hailey' or user_name == 'root':
        root = '/Users/hailey/Desktop/8_V_slab/figures'
    else:
        raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

    X_default = [
        'bulkox', 'bulk_e', 'bulk_ie', 'bulk_ieo', 'bulk_ien', 'bulk_cfse',
        'surfox', 'surf_e', 'surf_ie', 'surf_ieo', 'surf_ien', 'surf_cfse',
        'ICOHP', 'ICOHPo', 'ICOHPn', 'Hform', 'Hsub', 'Hevap', 'Hfus',
        'CN', 'group', 'period', 'chg', 'mag', 'volume', 'l_bond', 'cfse',
    ]
    ylabels = {
        'o': 'O Adsorption Energy (eV)',
        'oh': 'OH Adsorption Energy (eV)',
    }

    parser = argparse.ArgumentParser(description='Surface adsorption prediction using GPR, GBR, RF, LR, XGBoost, or LightGBM')
    parser.add_argument('--model', type=str, choices=['gpr', 'gbr', 'rf', 'lr', 'xgb', 'lgb'], default='gpr',
                      help='Model type to use (gpr, gbr, rf, lr, xgb, or lgb)')
    parser.add_argument('--Y', type=str, choices=['o', 'oh'], default='o',
                      help='Target column from lr_slab_data.csv (o_energy or oh_energy)')
    parser.add_argument('--X', nargs='+', default = X_default, help='List of feature columns from lr_slab_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size (default: 0.2)')
    parser.add_argument('--random_state', type=int, default=50, help='Random state for reproducibility (default: 50)')
    parser.add_argument('--energy_threshold', type=float, default=0.2,
                       help='Energy threshold (eV) for considering multiple coordinations in preference analysis')
    parser.add_argument('--corr_threshold', type=float, default=1.0,
                       help='Correlation threshold for feature selection (0.7, 0.8, 0.9 등)')
    parser.add_argument('--edge', action='store_true', help='only if n_electrons 0~10 (default: off)')
    parser.add_argument('--save', action='store_true', help='Save results to file')
    args = parser.parse_args()
    
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    YY = args.Y + '_energy'
    output_suffix = args.output

    row_str = ''.join(sorted(args.row)) if args.row else 'all'
    log_path = os.path.join(root, f'pred_{args.Y}_{args.model}_{row_str}_{output_suffix}.log')
    tsv_path = os.path.join(root, f'pred_{args.Y}_{args.model}_{row_str}_{output_suffix}.tsv')
    png_path = os.path.join(root, f'pred_{args.Y}_{args.model}_{row_str}_{output_suffix}.png')
    json_path = os.path.join(root, f'pred_{args.Y}_{args.model}_{row_str}_{output_suffix}.json')
    importance_png_path = os.path.join(root, f'pred_{args.Y}_{args.model}_{row_str}_{output_suffix}_importance.png')

    try:
        # Load data
        print("Loading data...")
        load_start = time.time()
        df_slab = pd.read_csv(os.path.join(root, 'lr_slab_data.csv'), index_col=0)
        print_time("Data loading completed", time.time() - load_start)
    except FileNotFoundError as e:
        print(f"{RED}Error: Required data file not found: {e}{ENDC}")
        exit(1)

    print("Preprocessing data...")
    preprocess_start = time.time()
    df = df_slab.copy()
    df = df.dropna(subset=args.X + [YY])

    if args.row:
        df = df[df['row'].isin(args.row)]
    if args.coord:
        df = df[df['coord'].isin(args.coord)]
    if args.edge:
        df = df[(df['n_electrons'] >= 0) & (df['n_electrons'] <= 10)]

    # Handle missing values in both X and y
    X = df[args.X].astype(float)
    y = df[YY].astype(float)

    # Handle remaining missing values in X
    imputer = SimpleImputer(strategy='median')
    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns, index=X.index)

    # Calculate correlation matrix
    print("Calculating correlation matrix...")
    corr_matrix = df[args.X + [YY]].corr()
     
    # Update X with imputed features
    X = X_imputed
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.random_state
    )

    # Print training and test dataset sizes
    print(f"\nTraining and test dataset sizes:")
    print(f"Training data: {len(X_train)} samples")
    print(f"Test data: {len(X_test)} samples")
    print(f"Total data: {len(X_train) + len(X_test)} samples")

    # Scale the features
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

    # Define and fit model
    print(f"Setting up {args.model.upper()} model...")
    model_start = time.time()
    
    if args.model == 'gpr':
        # GPR model setup
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

        # Bayesian Optimization for GPR
        search_space = {
            'kernel__k1__k2__length_scale': Real(1e-3, 1e5, prior='log-uniform'),
            'kernel__k1__k2__alpha': Real(1e-3, 1e5, prior='log-uniform'),
            'kernel__k2__noise_level': Real(1e-3, 1e5, prior='log-uniform')
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
    elif args.model == 'lr':
        # LR model setup
        model = LinearRegression()
        # No hyperparameter optimization needed for LR
        search_space = {}
    elif args.model == 'xgb':
        # XGBoost model setup
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
        search_space = {
            'n_estimators': Integer(100, 300),
            'learning_rate': Real(0.01, 0.1),
            'max_depth': Integer(3, 8),
            'min_child_weight': Integer(1, 5),
            'subsample': Real(0.6, 0.9),
            'colsample_bytree': Real(0.6, 0.9),
            'gamma': Real(0, 0.5),
            'reg_alpha': Real(0, 1),
            'reg_lambda': Real(0, 1)
        }
    else:  # lgb
        # LightGBM model setup
        model = lgb.LGBMRegressor(
            n_estimators=1000,
            learning_rate=0.01,
            num_leaves=31,
            random_state=args.random_state,
            n_jobs=1,  # 단일 프로세스 사용
            verbose=-1  # 경고 메시지 출력 억제
        )
        
        # Bayesian Optimization for LightGBM
        search_space = {
            'n_estimators': Integer(100, 2000),
            'learning_rate': Real(0.001, 0.1, prior='log-uniform'),
            'num_leaves': Integer(20, 100),
            'min_child_samples': Integer(5, 100),
            'subsample': Real(0.6, 1.0),
            'colsample_bytree': Real(0.6, 1.0)
        }
        
        # Ensure X_train_scaled and X_test_scaled are DataFrames with feature names
        X_train_scaled_values = X_train_scaled.values
        X_test_scaled_values = X_test_scaled.values
        
        bayes_search = BayesSearchCV(
            model,
            search_space,
            n_iter=50,
            cv=5,
            n_jobs=1,  # 단일 프로세스 사용
            random_state=args.random_state,
            verbose=0
        )
        
        # Fit the model with feature names
        bayes_search.fit(X_train_scaled, y_train)
        model = bayes_search.best_estimator_
        
        # Make predictions with feature names
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        std_train = get_prediction_std(model, X_train_scaled, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled, model_type=args.model)

    # Perform cross-validation
    print("Performing cross-validation...")
    cv_start = time.time()
    kf = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    
    # 경고 메시지를 방지하기 위해 NumPy 배열로 변환
    if args.model == 'lgb':
        # LightGBM의 경우 DataFrame 유지
        X_train_scaled_values = X_train_scaled
    else:
        X_train_scaled_values = X_train_scaled.values
    
    # Calculate multiple metrics for cross-validation
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
    
    # 모델 학습 추가
    print("Fitting model...")
    model.fit(X_train_scaled_values, y_train)
    
    # Make predictions
    print("Making predictions...")
    pred_start = time.time()
    
    # 예측을 위해 NumPy 배열로 변환
    X_train_scaled_values = X_train_scaled.values
    X_test_scaled_values = X_test_scaled.values
    
    if args.model == 'gpr':
        y_pred_train, std_train = model.predict(X_train_scaled_values, return_std=True)
        y_pred_test, std_test = model.predict(X_test_scaled_values, return_std=True)
    elif args.model == 'gbr':
        y_pred_train = model.predict(X_train_scaled_values)
        y_pred_test = model.predict(X_test_scaled_values)
        std_train = get_prediction_std(model, X_train_scaled_values, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled_values, model_type=args.model)
    elif args.model == 'rf':
        y_pred_train = model.predict(X_train_scaled_values)
        y_pred_test = model.predict(X_test_scaled_values)
        std_train = get_prediction_std(model, X_train_scaled_values, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled_values, model_type=args.model)
    elif args.model == 'lr':
        y_pred_train = model.predict(X_train_scaled_values)
        y_pred_test = model.predict(X_test_scaled_values)
        std_train = get_prediction_std(model, X_train_scaled_values, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled_values, model_type=args.model)
    elif args.model == 'lgb':
        # LightGBM의 경우 feature names를 유지
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        std_train = get_prediction_std(model, X_train_scaled, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled, model_type=args.model)
    else:  # xgb
        y_pred_train = model.predict(X_train_scaled_values)
        y_pred_test = model.predict(X_test_scaled_values)
        std_train = get_prediction_std(model, X_train_scaled_values, model_type=args.model)
        std_test = get_prediction_std(model, X_test_scaled_values, model_type=args.model)
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
    importance = feature_importance(X_train_scaled_values, y_train, model, args.X, args.model)
    cumulative_scores = plot_feature_importance(X_train_scaled, X_test_scaled, y_train, y_test, args.X, importance_png_path, args.model)
    print_time("Feature importance calculation completed", time.time() - importance_start)

    # Save results
    print("Saving results...")
    save_start = time.time()

    # Save metrics and cross-validation results
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
        
        # Add Feature Importance Analysis
        f.write("\nFeature Importance Analysis:\n")
        for score in cumulative_scores:
            f.write(f"\nNumber of features: {score['n_features']}\n")
            f.write(f"Features used: {', '.join(score['features'])}\n")
            f.write(f"Train MAE: {score['mae_train']:.4f}\n")
            f.write(f"Test MAE: {score['mae_test']:.4f}\n")
    print(f"{BLUE}Log file saved as {log_path}{ENDC}")

    # Save predictions
    print("Saving predictions...")
    
    # 전체 데이터에 대한 예측을 한 번에 수행
    X_all = X  # 이미 imputed된 데이터 사용
    X_all_scaled = pd.DataFrame(
        scaler.transform(X_all),
        columns=X_all.columns,
        index=X_all.index
    )
    
    # LightGBM 모델의 경우 feature names 유지
    if args.model == 'lgb':
        y_pred_all = model.predict(X_all_scaled)
        std_all = get_prediction_std(model, X_all_scaled, model_type=args.model)
    else:
        # 다른 모델의 경우 NumPy 배열로 변환
        X_all_scaled_values = X_all_scaled.values
        y_pred_all = model.predict(X_all_scaled_values)
        std_all = get_prediction_std(model, X_all_scaled_values, model_type=args.model)
    
    # 예측 결과를 DataFrame으로 저장
    predictions_df = pd.DataFrame({
        'Y_true': df[YY].values,
        'Y_pred': y_pred_all,
        'row': df['row'].values,
        'coord': df['coord'].values,
        'metal': df['metal'].values
    })
    
    # DataFrame 생성
    df_result = pd.DataFrame({
        'metal': df['metal'].values,
        'row': df['row'].values,
        'coord': df['coord'].values,
        'Y_true': y.values,
        'Y_pred': y_pred_all,
        'std': std_all
    })
    
    df_result.to_csv(tsv_path, sep='\t', index=False)
    print(f"{BLUE}TSV file saved as {tsv_path}{ENDC}")

    # Plot parity with color by 'row' and marker by 'coord'
    row_map = {'3d': 'red', '4d': 'green', '5d': 'blue'}
    coord_map = {'WZ': '+', 'ZB': 'x', 'TN': 'o', 'PD': 'o', 'NB': 's', 'RS': 'D', 'LT': 'h', '+3': 'v', '+4': '^', '+5': '<', '+6': '>'}

    plt.figure(figsize=(10, 6))
    
    # 이미 imputed된 데이터 사용
    X_all = X
    X_all_scaled = pd.DataFrame(
        scaler.transform(X_all),
        columns=X_all.columns,
        index=X_all.index
    )
    
    if args.model == 'lgb':
        y_pred_all = model.predict(X_all_scaled)
    else:
        X_all_scaled_values = X_all_scaled.values
        y_pred_all = model.predict(X_all_scaled_values)
    
    predictions_df = pd.DataFrame({
        'Y_true': df[YY].values,
        'Y_pred': y_pred_all,
        'row': df['row'],
        'coord': df['coord'],
        'metal': df['metal']
    })
    
    # 각 row와 coordination type에 대해 플롯
    for r in df['row'].unique():
        for c in df['coord'].unique():
            subset = predictions_df[
                (predictions_df['row'] == r) & 
                (predictions_df['coord'] == c)
            ]
            if subset.empty:
                continue
            
            plt.scatter(
                subset['Y_true'],
                subset['Y_pred'],
                label=f'{r}_{c}',
                alpha=0.3,
                color=row_map.get(r, 'gray'),
                marker=coord_map.get(c, 'x')
            )
            
            # 금속 이름 표시
            for _, row_data in subset.iterrows():
                plt.annotate(
                    row_data['metal'], 
                    (row_data['Y_true'], row_data['Y_pred']), 
                    fontsize=8
                )

    plt.plot([y.min(), y.max()], [y.min(), y.max()], '--', lw=1, color='black')
    plt.xlabel(f'DFT-calculated {ylabels[args.Y]}')
    plt.ylabel(f'Predicted {ylabels[args.Y]}')
    plt.legend(bbox_to_anchor=(-0.10, 1.01), loc='upper right', labelspacing=0.3, fontsize=9)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(png_path, transparent=True, dpi=300)
    plt.close()
    print(f"{BLUE}Parity plot saved as {png_path}{ENDC}")

    # Save all results in JSON format
    results = {
        'metrics': metrics,
        'model_params': {
            'type': args.model,
            'best_params': getattr(bayes_search, 'best_params_', None) if 'bayes_search' in locals() else None
        },
        'feature_names': args.X,
        'target': YY,
        'feature_importance': importance.to_dict()
    }
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"{BLUE}Results saved as {json_path}{ENDC}")

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
        if args.model in ['gpr', 'gbr', 'rf'] and 'bayes_start' in locals():
            f.write(f"Bayesian Optimization: {time.time() - bayes_start:.2f} seconds\n")
        f.write(f"Cross-validation: {time.time() - cv_start:.2f} seconds\n")
        f.write(f"Predictions: {time.time() - pred_start:.2f} seconds\n")
        f.write(f"Metrics calculation: {time.time() - metrics_start:.2f} seconds\n")
        f.write(f"Feature importance calculation: {time.time() - importance_start:.2f} seconds\n")
        f.write(f"Results saving and plotting: {time.time() - save_start:.2f} seconds\n")

    # Save results to file if --save option is used
    if args.save:
        result_filename = f'pred_{args.Y}_{args.model}_all_result.log'
        
        with open(result_filename, 'w') as f:
            f.write(f'Model: {args.model}\n')
            f.write(f'Features: {", ".join(args.X)}\n')  # args.X 사용
            f.write(f'Target: {YY}\n')
            f.write(f'R2 Score: {metrics["test"]["r2"]:.4f}\n')
            f.write(f'MAE: {metrics["test"]["mae"]:.4f}\n')
            f.write(f'RMSE: {metrics["test"]["rmse"]:.4f}\n')
            f.write(f'MAPE: {metrics["test"]["max_error"] / y.max() * 100:.4f}%\n')
            f.write(f'Max Error: {metrics["test"]["max_error"]:.4f}\n')
            f.write(f'Training Time: {time.time() - model_start:.2f} seconds\n')
            f.write(f'Prediction Time: {time.time() - pred_start:.2f} seconds\n')
            f.write('\nFeature Importance:\n')
            for feature, imp_value in importance.items():
                f.write(f'{feature}: {imp_value:.4f}\n')

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
    """Calculate feature importance (permutation importance, MAE 기준)"""
    if model_type == 'gpr':
        kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e5)) * RationalQuadratic(
            length_scale=1.0,
            alpha=1.0,
            length_scale_bounds=(1e-3, 1e6),
            alpha_bounds=(1e-3, 1e6)
        ) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-3, 1e5))
        
        # DataFrame을 numpy 배열로 변환
        if isinstance(X, pd.DataFrame):
            X_array = X.values
        else:
            X_array = X
        
        # 원래 MAE 계산
        base_mae = mean_absolute_error(y, model.predict(X_array))
        importance = []
        for i in tqdm(range(X_array.shape[1]), desc="Calculating feature importance (MAE)"):
            X_permuted = X_array.copy()
            X_permuted[:, i] = np.random.permutation(X_permuted[:, i])
            permuted_mae = mean_absolute_error(y, model.predict(X_permuted))
            importance.append(permuted_mae - base_mae)
        return pd.Series(importance, index=feature_names)
    elif model_type in ['gbr', 'rf']:
        if isinstance(X, pd.DataFrame):
            X_array = X.values
        else:
            X_array = X
        
        from sklearn.base import clone
        temp_model = clone(model)
        temp_model.fit(X_array, y)
        
        importance = temp_model.feature_importances_
        return pd.Series(importance, index=feature_names)
    elif model_type == 'xgb':
        if isinstance(X, pd.DataFrame):
            X_array = X.values
        else:
            X_array = X
        
        temp_model = xgb.XGBRegressor()
        temp_model.fit(X_array, y)
        importance = temp_model.feature_importances_
        return pd.Series(importance, index=feature_names)
    elif model_type == 'lgb':
        if isinstance(X, pd.DataFrame):
            X_array = X
        else:
            X_array = pd.DataFrame(X, columns=feature_names)
        
        temp_model = lgb.LGBMRegressor()
        temp_model.fit(X_array, y)
        importance = temp_model.feature_importances_
        return pd.Series(importance, index=feature_names)
    else:  # lr
        if isinstance(X, pd.DataFrame):
            X_array = X.values
        else:
            X_array = X
        
        importance = np.abs(model.coef_)
        return pd.Series(importance, index=feature_names)

def get_prediction_std(model, X, n_iterations=100, model_type='gpr'):
    """
    Estimate prediction uncertainty
    """
    # X를 numpy 배열로 변환하되, feature names 유지
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
        X_array = X
    else:
        if hasattr(model, 'feature_names_in_'):
            feature_names = model.feature_names_in_
            X_array = pd.DataFrame(X, columns=feature_names)
        else:
            X_array = X

    if model_type == 'gpr':
        _, std = model.predict(X_array, return_std=True)
        return std
    elif model_type in ['gbr', 'rf']:
        predictions = []
        n_estimators = len(model.estimators_) if hasattr(model, 'estimators_') else model.n_estimators
        
        for i in range(n_iterations):
            if model_type == 'rf':
                indices = np.random.choice(n_estimators, size=n_estimators//2, replace=False)
                pred = np.zeros(X_array.shape[0])
                for idx in indices:
                    pred += model.estimators_[idx].predict(X_array)
                pred /= len(indices)
            else:  # GBR
                pred = np.zeros(X_array.shape[0])
                for y_pred in model.staged_predict(X_array):
                    pred = y_pred
                predictions.append(pred)
                continue
            predictions.append(pred)
        return np.std(predictions, axis=0)
    elif model_type == 'lgb':
        predictions = []
        for i in range(n_iterations):
            # X_array가 DataFrame인지 확인하고, 아니면 DataFrame으로 변환
            if not isinstance(X_array, pd.DataFrame):
                X_array = pd.DataFrame(X_array, columns=model.feature_names_in_)
            pred = model.predict(X_array)
            predictions.append(pred)
        return np.std(predictions, axis=0)
    else:  # lr
        y_pred = model.predict(X_array)
        residuals = y_pred - model.predict(X_array)
        return np.std(residuals) * np.ones(X_array.shape[0])

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
        X_train_selected = X_train[[feature]]
        X_test_selected = X_test[[feature]]
        
        # Scale the selected features
        scaler_selected = RobustScaler()
        X_train_scaled_selected = pd.DataFrame(
            scaler_selected.fit_transform(X_train_selected),
            columns=X_train_selected.columns
        )
        X_test_scaled_selected = pd.DataFrame(
            scaler_selected.transform(X_test_selected),
            columns=X_test_selected.columns
        )
        
        # Train model with selected features
        if model_type == 'gpr':
            kernel = ConstantKernel(1.0) * RationalQuadratic() + WhiteKernel()
            model_selected = MyGPR(kernel=kernel, random_state=42)
        elif model_type == 'lr':
            model_selected = LinearRegression()
        else:  # gbr
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
    
    # Sort features by individual MAE
    sorted_features = sorted(single_feature_scores, key=lambda x: x['mae_train'])
    sorted_feature_names = [score['feature'] for score in sorted_features]
    
    # Calculate cumulative MAE
    for i in range(1, len(sorted_feature_names) + 1):
        selected_features = sorted_feature_names[:i]
        X_train_selected = X_train[selected_features]
        X_test_selected = X_test[selected_features]
        
        # Scale the selected features
        scaler_selected = RobustScaler()
        X_train_scaled_selected = pd.DataFrame(
            scaler_selected.fit_transform(X_train_selected),
            columns=selected_features
        )
        X_test_scaled_selected = pd.DataFrame(
            scaler_selected.transform(X_test_selected),
            columns=selected_features
        )
        
        # Train model with selected features
        if model_type == 'gpr':
            kernel = ConstantKernel(1.0) * RationalQuadratic() + WhiteKernel()
            model_selected = MyGPR(kernel=kernel, random_state=42)
        elif model_type == 'lr':
            model_selected = LinearRegression()
        else:  # gbr
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

if __name__ == '__main__':
    print("Script started...")
    print("Importing required libraries...")
    main() 