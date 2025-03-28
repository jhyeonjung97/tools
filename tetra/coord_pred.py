import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from tqdm import tqdm
import time
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.pipeline import Pipeline
from skopt import BayesSearchCV
from skopt.space import Real, Integer
import socket
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
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

def feature_importance(X, y, model, feature_names):
    """Calculate feature importance for classification model"""
    if hasattr(model, 'feature_importances_'):
        importance = pd.Series(model.feature_importances_, index=feature_names)
    else:
        importance = pd.Series(0, index=feature_names)
    return importance

def plot_confusion_matrix(y_true, y_pred, output_path):
    """Plot confusion matrix"""
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_feature_importance(importance, output_path):
    """Plot feature importance"""
    plt.figure(figsize=(12, 6))
    importance.sort_values(ascending=True).plot(kind='barh')
    plt.xlabel('Feature Importance')
    plt.title('Feature Importance Analysis')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def analyze_coordination_types(df_result):
    """Analyze coordination type predictions"""
    def get_coordination_type(coord):
        if coord in ['WZ', 'ZB']:
            return 'tetrahedral'
        elif coord in ['TN', 'PD', 'NB']:
            return 'squareplanar'
        elif coord == 'RS':
            return 'octahedral'
        elif coord == 'LT':
            return 'pyramidal'
        return 'unknown'
    
    df_result['true_type'] = df_result['coord_true'].apply(get_coordination_type)
    df_result['pred_type'] = df_result['coord_pred'].apply(get_coordination_type)
    df_result['type_match'] = df_result['true_type'] == df_result['pred_type']
    
    return df_result

def main():
    start_time = time.time()
    print("Starting coordination prediction analysis...")
    
    parser = argparse.ArgumentParser(description='Coordination prediction using bulk_data.csv and mendeleev_data.csv')
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
    parser.add_argument('--model', type=str, default='rf', choices=['rf', 'gb', 'lr', 'svm', 'xgb', 'lgb'], help='Model type: Random Forest (rf), Gradient Boosting (gb), Logistic Regression (lr), Support Vector Machine (svm), XGBoost (xgb), LightGBM (lgb)')
    args = parser.parse_args()
    
    # Convert feature names if they start with ICOHP or ICOOP (prepend '-')
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    
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
    df = df.dropna(subset=args.X + ['coord'])

    X = df[args.X].astype(float)
    y = df['coord']

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=args.test_size, random_state=args.random_state)

    # Scale the features
    print("Scaling features...")
    scale_start = time.time()
    scaler = RobustScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    print_time("Feature scaling completed", time.time() - scale_start)

    # Define and fit classification model
    print("Setting up classification model...")
    model_start = time.time()
    
    if args.model == 'rf':
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=None,
            min_samples_split=2,
            min_samples_leaf=1,
            random_state=args.random_state
        )
    elif args.model == 'gb':
        model = GradientBoostingClassifier(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=3,
            random_state=args.random_state
        )
    elif args.model == 'lr':
        model = LogisticRegression(max_iter=1000, random_state=args.random_state)
    elif args.model == 'svm':
        model = SVC(random_state=args.random_state, probability=True)
    elif args.model == 'xgb':
        model = xgb.XGBClassifier(n_estimators=100, random_state=args.random_state)
    else:  # gb
        model = GradientBoostingClassifier(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=3,
            random_state=args.random_state
        )
    elif args.model == 'lgb':
        model = lgb.LGBMClassifier(n_estimators=100, random_state=args.random_state)

    # Perform Bayesian Optimization
    print("Performing Bayesian Optimization...")
    grid_start = time.time()
    if args.model == 'rf':
        search_space = {
            'n_estimators': Integer(50, 200),
            'max_depth': Integer(5, 30),
            'min_samples_split': Integer(2, 10),
            'min_samples_leaf': Integer(1, 5)
        }
    elif args.model == 'gb':
        search_space = {
            'n_estimators': Integer(50, 200),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'max_depth': Integer(3, 10)
        }
    elif args.model == 'lr':
        search_space = {
            'C': Real(1e-3, 1e3, prior='log-uniform'),
            'l1_ratio': Real(0, 1),
            'class_weight': ['balanced', None]
        }
    elif args.model == 'svm':
        search_space = {
            'C': Real(1e-3, 1e3, prior='log-uniform'),
            'gamma': Real(1e-4, 1e-1, prior='log-uniform'),
            'kernel': ['rbf', 'linear']
        }
    elif args.model == 'xgb':
        search_space = {
            'n_estimators': Integer(50, 200),
            'max_depth': Integer(3, 10),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'subsample': Real(0.6, 1.0),
            'colsample_bytree': Real(0.6, 1.0)
        }
    elif args.model == 'lgb':
        search_space = {
            'n_estimators': Integer(50, 200),
            'max_depth': Integer(3, 10),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'num_leaves': Integer(20, 100),
            'feature_fraction': Real(0.6, 1.0),
            'bagging_fraction': Real(0.6, 1.0)
        }
    else:  # gb
        search_space = {
            'n_estimators': Integer(50, 200),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'max_depth': Integer(3, 10)
        }

    bayes_search = BayesSearchCV(
        model,
        search_space,
        n_iter=50,
        cv=5,
        scoring='accuracy',
        n_jobs=1,
        random_state=args.random_state,
        verbose=1,
        n_points=5
    )
    bayes_search.fit(X_train_scaled, y_train)
    print_time("Bayesian Optimization completed", time.time() - grid_start)
    print(f"{MAGENTA}Best parameters: {bayes_search.best_params_}{ENDC}")
    
    # Use best model
    model = bayes_search.best_estimator_
    
    # Perform cross-validation
    print("Performing cross-validation...")
    cv_start = time.time()
    kf = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    
    # Calculate multiple metrics for cross-validation
    cv_accuracy = cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='accuracy')
    cv_precision = cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='precision_macro')
    cv_recall = cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='recall_macro')
    cv_f1 = cross_val_score(model, X_train_scaled, y_train, cv=kf, scoring='f1_macro')
    
    print_time("Cross-validation completed", time.time() - cv_start)
    print(f"{MAGENTA}Cross-validation Accuracy scores: {cv_accuracy}{ENDC}")
    print(f"{MAGENTA}Mean CV Accuracy score: {cv_accuracy.mean():.4f} (+/- {cv_accuracy.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation Precision scores: {cv_precision}{ENDC}")
    print(f"{MAGENTA}Mean CV Precision score: {cv_precision.mean():.4f} (+/- {cv_precision.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation Recall scores: {cv_recall}{ENDC}")
    print(f"{MAGENTA}Mean CV Recall score: {cv_recall.mean():.4f} (+/- {cv_recall.std() * 2:.4f}){ENDC}")
    print(f"{MAGENTA}Cross-validation F1 scores: {cv_f1}{ENDC}")
    print(f"{MAGENTA}Mean CV F1 score: {cv_f1.mean():.4f} (+/- {cv_f1.std() * 2:.4f}){ENDC}")
    
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
            'accuracy': accuracy_score(y_train, y_pred_train),
            'precision': precision_score(y_train, y_pred_train, average='macro'),
            'recall': recall_score(y_train, y_pred_train, average='macro'),
            'f1': f1_score(y_train, y_pred_train, average='macro')
        },
        'test': {
            'accuracy': accuracy_score(y_test, y_pred_test),
            'precision': precision_score(y_test, y_pred_test, average='macro'),
            'recall': recall_score(y_test, y_pred_test, average='macro'),
            'f1': f1_score(y_test, y_pred_test, average='macro')
        },
        'cv': {
            'accuracy': {
                'mean': cv_accuracy.mean(),
                'std': cv_accuracy.std(),
                'scores': cv_accuracy.tolist()
            },
            'precision': {
                'mean': cv_precision.mean(),
                'std': cv_precision.std(),
                'scores': cv_precision.tolist()
            },
            'recall': {
                'mean': cv_recall.mean(),
                'std': cv_recall.std(),
                'scores': cv_recall.tolist()
            },
            'f1': {
                'mean': cv_f1.mean(),
                'std': cv_f1.std(),
                'scores': cv_f1.tolist()
            }
        }
    }
    print_time("Metrics calculation completed", time.time() - metrics_start)

    # Calculate feature importance
    print("Calculating feature importance...")
    importance_start = time.time()
    importance = feature_importance(X_train_scaled, y_train, model, args.X)
    print_time("Feature importance calculation completed", time.time() - importance_start)

    # Save results
    output_suffix = args.output
    print("Saving results...")
    save_start = time.time()

    # Define paths for all output files
    log_path = os.path.join(root, f'coord_pred_{output_suffix}.log')
    tsv_path = os.path.join(root, f'coord_pred_{output_suffix}.tsv')
    json_path = os.path.join(root, f'coord_pred_{output_suffix}.json')
    cm_path = os.path.join(root, f'coord_pred_{output_suffix}_cm.png')
    importance_path = os.path.join(root, f'coord_pred_{output_suffix}_importance.png')

    # Save metrics and cross-validation results
    with open(log_path, 'w') as f:
        f.write("Training Metrics:\n")
        f.write(f"Accuracy: {metrics['train']['accuracy']:.4f}\n")
        f.write(f"Precision: {metrics['train']['precision']:.4f}\n")
        f.write(f"Recall: {metrics['train']['recall']:.4f}\n")
        f.write(f"F1 Score: {metrics['train']['f1']:.4f}\n\n")
        f.write("Test Metrics:\n")
        f.write(f"Accuracy: {metrics['test']['accuracy']:.4f}\n")
        f.write(f"Precision: {metrics['test']['precision']:.4f}\n")
        f.write(f"Recall: {metrics['test']['recall']:.4f}\n")
        f.write(f"F1 Score: {metrics['test']['f1']:.4f}\n\n")
        f.write("Cross-validation Results:\n")
        f.write(f"Accuracy scores: {cv_accuracy}\n")
        f.write(f"Mean Accuracy score: {cv_accuracy.mean():.4f} (+/- {cv_accuracy.std() * 2:.4f})\n")
        f.write(f"Precision scores: {cv_precision}\n")
        f.write(f"Mean Precision score: {cv_precision.mean():.4f} (+/- {cv_precision.std() * 2:.4f})\n")
        f.write(f"Recall scores: {cv_recall}\n")
        f.write(f"Mean Recall score: {cv_recall.mean():.4f} (+/- {cv_recall.std() * 2:.4f})\n")
        f.write(f"F1 scores: {cv_f1}\n")
        f.write(f"Mean F1 score: {cv_f1.mean():.4f} (+/- {cv_f1.std() * 2:.4f})\n")
    print(f"{BLUE}Log file saved as {log_path}{ENDC}")

    # Save predictions
    df_result = pd.DataFrame({
        'metal': df.index,
        'row': df['row'],
        'coord_true': y,
        'coord_pred': model.predict(scaler.transform(X))
    })
    
    # Extract actual metal name from the index
    df_result['metal'] = df_result.apply(lambda row: df.loc[row['metal'], 'metal'], axis=1)
    
    # Analyze coordination types
    df_result = analyze_coordination_types(df_result)
    
    # Add coordination type analysis to log file
    with open(log_path, 'a') as f:
        f.write("\nCoordination Type Analysis:\n")
        f.write("-------------------------\n")
        f.write(f"Overall type prediction accuracy: {df_result['type_match'].mean():.4f}\n\n")
        
        f.write("Accuracy by coordination type:\n")
        for coord_type in ['tetrahedral', 'squareplanar', 'octahedral', 'pyramidal']:
            type_data = df_result[df_result['true_type'] == coord_type]
            if len(type_data) > 0:
                f.write(f"\n{coord_type.capitalize()}:\n")
                f.write(f"Number of samples: {len(type_data)}\n")
                f.write(f"Accuracy: {type_data['type_match'].mean():.4f}\n")
        
        f.write("\nAccuracy by row:\n")
        for row in ['3d', '4d', '5d']:
            row_data = df_result[df_result['row'] == row]
            if len(row_data) > 0:
                f.write(f"\n{row} metals:\n")
                f.write(f"Number of samples: {len(row_data)}\n")
                f.write(f"Accuracy: {row_data['type_match'].mean():.4f}\n")
    
    df_result.to_csv(tsv_path, sep='\t', index=False)
    print(f"{BLUE}TSV file saved as {tsv_path}{ENDC}")

    # Plot confusion matrix
    plot_confusion_matrix(y_test, y_pred_test, cm_path)
    print(f"{BLUE}Confusion matrix plot saved as {cm_path}{ENDC}")

    # Plot feature importance
    plot_feature_importance(importance, importance_path)
    print(f"{BLUE}Feature importance plot saved as {importance_path}{ENDC}")

    # Save all results in JSON format
    results = {
        'metrics': metrics,
        'model_params': {
            'model_type': args.model,
            'best_params': bayes_search.best_params_,
            'random_state': args.random_state
        },
        'feature_names': args.X,
        'feature_importance': importance.to_dict()
    }
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"{BLUE}Results saved as {json_path}{ENDC}")

    # Print total execution time
    total_time = time.time() - start_time
    print(f"Execution Summary:")
    print_time("Total execution time", total_time)
    print_time("Model training and evaluation time", time.time() - model_start)
    
    # Save execution time to log file
    with open(log_path, 'a') as f:
        f.write("\nExecution Times:\n")
        f.write(f"Total execution time: {total_time:.2f} seconds\n")
        f.write(f"Data loading: {time.time() - load_start:.2f} seconds\n")
        f.write(f"Data preprocessing: {time.time() - preprocess_start:.2f} seconds\n")
        f.write(f"Feature scaling: {time.time() - scale_start:.2f} seconds\n")
        f.write(f"Bayesian Optimization: {time.time() - grid_start:.2f} seconds\n")
        f.write(f"Cross-validation: {time.time() - cv_start:.2f} seconds\n")
        f.write(f"Predictions: {time.time() - pred_start:.2f} seconds\n")
        f.write(f"Metrics calculation: {time.time() - metrics_start:.2f} seconds\n")
        f.write(f"Feature importance calculation: {time.time() - importance_start:.2f} seconds\n")
        f.write(f"Results saving and plotting: {time.time() - save_start:.2f} seconds\n")

if __name__ == '__main__':
    main() 