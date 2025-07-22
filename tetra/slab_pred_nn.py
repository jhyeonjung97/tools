import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from tqdm import tqdm
import time
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, max_error
import socket
import warnings
warnings.filterwarnings('ignore')

# Neural Network imports
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import KFold

# Hyperparameter optimization
import optuna
from optuna.samplers import TPESampler

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

class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_sizes=[128, 64, 32], dropout_rate=0.2):
        super(NeuralNetwork, self).__init__()
        layers = []
        prev_size = input_size
        
        for hidden_size in hidden_sizes:
            layers.extend([
                nn.Linear(prev_size, hidden_size),
                nn.ReLU(),
                nn.Dropout(dropout_rate),
                nn.LayerNorm(hidden_size)  # safer for small batch sizes
            ])
            prev_size = hidden_size
        
        layers.append(nn.Linear(prev_size, 1))
        
        self.network = nn.Sequential(*layers)
        
    def forward(self, x):
        return self.network(x)

def train_neural_network(model, train_loader, val_loader, epochs=100, lr=0.001, patience=20, weight_decay=1e-5):
    """Train neural network with early stopping"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)
    
    best_val_loss = float('inf')
    patience_counter = 0
    train_losses = []
    val_losses = []
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0.0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_X).squeeze()
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X).squeeze()
                loss = criterion(outputs, batch_y)
                val_loss += loss.item()
        
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        train_losses.append(train_loss)
        val_losses.append(val_loss)
        
        scheduler.step(val_loss)
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            best_model_state = model.state_dict().copy()
        else:
            patience_counter += 1
            
        if patience_counter >= patience:
            print(f"Early stopping at epoch {epoch}")
            break
    
    # Load best model
    model.load_state_dict(best_model_state)
    return model, train_losses, val_losses

def predict_with_uncertainty(model, X, n_iterations=100):
    """Predict with uncertainty using dropout"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model.train()  # Enable dropout for uncertainty estimation
    
    predictions = []
    X_tensor = torch.FloatTensor(X).to(device)
    
    with torch.no_grad():
        for _ in range(n_iterations):
            pred = model(X_tensor).cpu().numpy().flatten()
            predictions.append(pred)
    
    predictions = np.array(predictions)
    mean_pred = np.mean(predictions, axis=0)
    std_pred = np.std(predictions, axis=0)
    
    return mean_pred, std_pred

def feature_importance_nn(model, X, y, feature_names, n_iterations=50):
    """Calculate feature importance using permutation importance"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model.eval()
    
    # Ensure X and y are numpy arrays
    X_np = X.values if isinstance(X, pd.DataFrame) else X
    y_np = y.values if isinstance(y, pd.Series) else y
    X_tensor = torch.FloatTensor(X_np).to(device)
    y_tensor = torch.FloatTensor(y_np).to(device)
    
    # Calculate baseline loss
    with torch.no_grad():
        baseline_pred = model(X_tensor).squeeze()
        baseline_loss = nn.MSELoss()(baseline_pred, y_tensor).item()
    
    importance = []
    for i in tqdm(range(X_tensor.shape[1]), desc="Calculating feature importance"):
        permuted_losses = []
        
        for _ in range(n_iterations):
            X_permuted = X_tensor.clone()
            X_permuted[:, i] = X_permuted[:, i][torch.randperm(X_permuted.size(0))]
            
            with torch.no_grad():
                permuted_pred = model(X_permuted).squeeze()
                permuted_loss = nn.MSELoss()(permuted_pred, y_tensor).item()
                permuted_losses.append(permuted_loss)
        
        importance.append(np.mean(permuted_losses) - baseline_loss)
    
    return pd.Series(importance, index=feature_names)

def plot_feature_importance_nn(X_train, X_test, y_train, y_test, feature_names, output_path, hidden_sizes, dropout_rate, learning_rate, batch_size, weight_decay):
    """Create a cumulative MAE plot for feature importance"""
    single_feature_scores = []
    cumulative_scores = []

    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train.values)
    X_test_tensor = torch.FloatTensor(X_test.values)
    y_train_tensor = torch.FloatTensor(y_train.values)
    y_test_tensor = torch.FloatTensor(y_test.values)

    for feature in feature_names:
        # Select single feature
        feature_idx = list(feature_names).index(feature)
        X_train_selected = X_train_tensor[:, feature_idx:feature_idx+1]
        X_test_selected = X_test_tensor[:, feature_idx:feature_idx+1]
        
        # Scale the selected features
        scaler_selected = StandardScaler()
        X_train_scaled = torch.FloatTensor(scaler_selected.fit_transform(X_train_selected))
        X_test_scaled = torch.FloatTensor(scaler_selected.transform(X_test_selected))
        
        # Train model with selected features
        model_selected = NeuralNetwork(input_size=1, hidden_sizes=[32, 16], dropout_rate=dropout_rate)
        
        # Create data loaders
        train_dataset = TensorDataset(X_train_scaled, y_train_tensor)
        test_dataset = TensorDataset(X_test_scaled, y_test_tensor)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        # Train model
        model_selected, _, _ = train_neural_network(model_selected, train_loader, test_loader, epochs=50, lr=learning_rate, weight_decay=weight_decay)
        
        # Make predictions
        model_selected.eval()
        with torch.no_grad():
            y_pred_train = model_selected(X_train_scaled).squeeze().numpy()
            y_pred_test = model_selected(X_test_scaled).squeeze().numpy()
        
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
        feature_indices = [list(feature_names).index(feat) for feat in selected_features]
        
        X_train_selected = X_train_tensor[:, feature_indices]
        X_test_selected = X_test_tensor[:, feature_indices]
        
        # Scale the selected features
        scaler_selected = StandardScaler()
        X_train_scaled = torch.FloatTensor(scaler_selected.fit_transform(X_train_selected))
        X_test_scaled = torch.FloatTensor(scaler_selected.transform(X_test_selected))
        
        # Train model with selected features
        model_selected = NeuralNetwork(input_size=len(selected_features), hidden_sizes=[64, 32], dropout_rate=dropout_rate)
        
        # Create data loaders
        train_dataset = TensorDataset(X_train_scaled, y_train_tensor)
        test_dataset = TensorDataset(X_test_scaled, y_test_tensor)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        # Train model
        model_selected, _, _ = train_neural_network(model_selected, train_loader, test_loader, epochs=50, lr=learning_rate, weight_decay=weight_decay)
        
        # Make predictions
        model_selected.eval()
        with torch.no_grad():
            y_pred_train = model_selected(X_train_scaled).squeeze().numpy()
            y_pred_test = model_selected(X_test_scaled).squeeze().numpy()
        
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

def print_time(message, time_value):
    """Print time-related information with color coding"""
    print(f"{BOLD}{message}: {time_value:.2f} seconds{ENDC}")

def objective(trial, X_train, y_train, X_val, y_val, input_size):
    """Objective function for Optuna optimization"""
    # Suggest hyperparameters
    n_layers = trial.suggest_int('n_layers', 2, 4)
    hidden_sizes = []
    for i in range(n_layers):
        if i == 0:
            hidden_sizes.append(trial.suggest_int(f'hidden_size_{i}', 32, 256))
        else:
            hidden_sizes.append(trial.suggest_int(f'hidden_size_{i}', 16, hidden_sizes[i-1]))
    
    dropout_rate = trial.suggest_float('dropout_rate', 0.1, 0.5)
    learning_rate = trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True)
    batch_size = trial.suggest_categorical('batch_size', [16, 32, 64, 128])
    weight_decay = trial.suggest_float('weight_decay', 1e-6, 1e-3, log=True)
    
    # Create model
    model = NeuralNetwork(
        input_size=input_size,
        hidden_sizes=hidden_sizes,
        dropout_rate=dropout_rate
    )
    
    # Create data loaders
    train_dataset = TensorDataset(
        torch.FloatTensor(X_train.values),
        torch.FloatTensor(y_train.values)
    )
    val_dataset = TensorDataset(
        torch.FloatTensor(X_val.values),
        torch.FloatTensor(y_val.values)
    )
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    
    # Train model with early stopping
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)
    
    best_val_loss = float('inf')
    patience_counter = 0
    patience = 15  # Shorter patience for optimization
    
    for epoch in range(50):  # Fewer epochs for optimization
        # Training
        model.train()
        train_loss = 0.0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_X).squeeze()
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X).squeeze()
                loss = criterion(outputs, batch_y)
                val_loss += loss.item()
        
        val_loss /= len(val_loader)
        scheduler.step(val_loss)
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            
        if patience_counter >= patience:
            break
    
    return best_val_loss

def optimize_hyperparameters(X_train, y_train, input_size, n_trials=50):
    """Optimize hyperparameters using Optuna"""
    print("Starting hyperparameter optimization...")
    opt_start = time.time()
    
    # Split training data for optimization
    X_opt_train, X_opt_val, y_opt_train, y_opt_val = train_test_split(
        X_train, y_train, test_size=0.2, random_state=42
    )
    
    # Create study
    study = optuna.create_study(
        direction='minimize',
        sampler=TPESampler(seed=42)
    )
    
    # Optimize
    study.optimize(
        lambda trial: objective(trial, X_opt_train, y_opt_train, X_opt_val, y_opt_val, input_size),
        n_trials=n_trials
    )
    
    print_time("Hyperparameter optimization completed", time.time() - opt_start)
    print(f"{GREEN}Best trial: {study.best_trial.value:.6f}{ENDC}")
    print(f"{GREEN}Best parameters: {study.best_trial.params}{ENDC}")
    
    return study.best_trial.params

def main():
    start_time = time.time()
    print("Starting surface adsorption prediction analysis with Neural Network...")
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

    parser = argparse.ArgumentParser(description='Surface adsorption prediction using Neural Network')
    parser.add_argument('--Y', type=str, choices=['o', 'oh'], default='o',
                      help='Target column from lr_slab_data.csv (o_energy or oh_energy)')
    parser.add_argument('--X', nargs='+', default = X_default, help='List of feature columns from lr_slab_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    parser.add_argument('--test_size', type=float, default=0.2, help='Test set size (default: 0.2)')
    parser.add_argument('--random_state', type=int, default=50, help='Random state for reproducibility (default: 50)')
    parser.add_argument('--hidden_sizes', nargs='+', type=int, default=[128, 64, 32], 
                      help='Hidden layer sizes (default: 128 64 32)')
    parser.add_argument('--dropout_rate', type=float, default=0.2, help='Dropout rate (default: 0.2)')
    parser.add_argument('--epochs', type=int, default=100, help='Number of training epochs (default: 100)')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate (default: 0.001)')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size (default: 32)')
    parser.add_argument('--patience', type=int, default=20, help='Early stopping patience (default: 20)')
    parser.add_argument('--edge', action='store_true', help='n_electrons가 0~10인 데이터만 사용 (기본값: off)')
    parser.add_argument('--save', action='store_true', help='Save results to file')
    parser.add_argument('--no_optimize', action='store_true', help='Disable hyperparameter optimization (default: optimization enabled)')
    parser.add_argument('--n_trials', type=int, default=50, help='Number of optimization trials (default: 50)')
    args = parser.parse_args()
    
    args.X = [('-' + x if x.startswith('ICOHP') or x.startswith('ICOOP') else x) for x in args.X]
    YY = args.Y + '_energy'
    output_suffix = args.output

    row_str = ''.join(sorted(args.row)) if args.row else 'all'
    log_path = os.path.join(root, f'{args.Y}_pred_nn_{row_str}_{output_suffix}.log')
    tsv_path = os.path.join(root, f'{args.Y}_pred_nn_{row_str}_{output_suffix}.tsv')
    png_path = os.path.join(root, f'{args.Y}_pred_nn_{row_str}_{output_suffix}.png')
    json_path = os.path.join(root, f'{args.Y}_pred_nn_{row_str}_{output_suffix}.json')
    importance_png_path = os.path.join(root, f'{args.Y}_pred_nn_{row_str}_{output_suffix}_importance.png')

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
    # n_electrons 필터 옵션 적용
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

    # Define and train neural network
    print("Setting up Neural Network model...")
    model_start = time.time()
    
    # Optimize hyperparameters if requested
    if not args.no_optimize:
        best_params = optimize_hyperparameters(X_train_scaled, y_train, len(args.X), args.n_trials)
        
        # Extract optimized parameters
        n_layers = best_params['n_layers']
        hidden_sizes = [best_params[f'hidden_size_{i}'] for i in range(n_layers)]
        dropout_rate = best_params['dropout_rate']
        learning_rate = best_params['learning_rate']
        batch_size = best_params['batch_size']
        weight_decay = best_params['weight_decay']
        
        print(f"{CYAN}Using optimized hyperparameters:{ENDC}")
        print(f"  Hidden sizes: {hidden_sizes}")
        print(f"  Dropout rate: {dropout_rate:.3f}")
        print(f"  Learning rate: {learning_rate:.6f}")
        print(f"  Batch size: {batch_size}")
        print(f"  Weight decay: {weight_decay:.6f}")
    else:
        # Use default parameters
        hidden_sizes = args.hidden_sizes
        dropout_rate = args.dropout_rate
        learning_rate = args.lr
        batch_size = args.batch_size
        weight_decay = 1e-5
    
    # Create neural network
    model = NeuralNetwork(
        input_size=len(args.X),
        hidden_sizes=hidden_sizes,
        dropout_rate=dropout_rate
    )
    
    # Create data loaders
    train_dataset = TensorDataset(
        torch.FloatTensor(X_train_scaled.values),
        torch.FloatTensor(y_train.values)
    )
    test_dataset = TensorDataset(
        torch.FloatTensor(X_test_scaled.values),
        torch.FloatTensor(y_test.values)
    )
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    # Train model
    print("Training neural network...")
    model, train_losses, val_losses = train_neural_network(
        model, train_loader, test_loader, 
        epochs=args.epochs, lr=learning_rate, patience=args.patience, weight_decay=weight_decay
    )
    print_time("Model training completed", time.time() - model_start)

    # Perform cross-validation
    print("Performing cross-validation...")
    cv_start = time.time()
    kf = KFold(n_splits=5, shuffle=True, random_state=args.random_state)
    
    cv_r2_scores = []
    cv_mae_scores = []
    cv_mse_scores = []
    
    for train_idx, val_idx in kf.split(X_train_scaled):
        X_cv_train, X_cv_val = X_train_scaled.iloc[train_idx], X_train_scaled.iloc[val_idx]
        y_cv_train, y_cv_val = y_train.iloc[train_idx], y_train.iloc[val_idx]
        
        # Create data loaders for CV
        cv_train_dataset = TensorDataset(
            torch.FloatTensor(X_cv_train.values),
            torch.FloatTensor(y_cv_train.values)
        )
        cv_val_dataset = TensorDataset(
            torch.FloatTensor(X_cv_val.values),
            torch.FloatTensor(y_cv_val.values)
        )
        
        cv_train_loader = DataLoader(cv_train_dataset, batch_size=batch_size, shuffle=True)
        cv_val_loader = DataLoader(cv_val_dataset, batch_size=batch_size, shuffle=False)
        
        # Train model for this fold
        cv_model = NeuralNetwork(
            input_size=len(args.X),
            hidden_sizes=hidden_sizes,
            dropout_rate=dropout_rate
        )
        cv_model, _, _ = train_neural_network(
            cv_model, cv_train_loader, cv_val_loader,
            epochs=args.epochs, lr=learning_rate, patience=args.patience, weight_decay=weight_decay
        )
        
        # Make predictions
        cv_model.eval()
        with torch.no_grad():
            y_pred_cv = cv_model(torch.FloatTensor(X_cv_val.values)).squeeze().numpy()
        
        # Calculate metrics
        cv_r2_scores.append(r2_score(y_cv_val, y_pred_cv))
        cv_mae_scores.append(mean_absolute_error(y_cv_val, y_pred_cv))
        cv_mse_scores.append(mean_squared_error(y_cv_val, y_pred_cv))
    
    cv_r2 = np.array(cv_r2_scores)
    cv_mae = np.array(cv_mae_scores)
    cv_mse = np.array(cv_mse_scores)
    
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
    
    model.eval()
    with torch.no_grad():
        y_pred_train = model(torch.FloatTensor(X_train_scaled.values)).squeeze().numpy()
        y_pred_test = model(torch.FloatTensor(X_test_scaled.values)).squeeze().numpy()
    
    # Get uncertainty estimates
    y_pred_train_mean, std_train = predict_with_uncertainty(model, X_train_scaled.values)
    y_pred_test_mean, std_test = predict_with_uncertainty(model, X_test_scaled.values)
    
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
    importance = feature_importance_nn(model, X_train_scaled, y_train, args.X, 50)
    cumulative_scores = plot_feature_importance_nn(X_train_scaled, X_test_scaled, y_train, y_test, args.X, importance_png_path, hidden_sizes, dropout_rate, learning_rate, batch_size, weight_decay)
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
    
    # Make predictions with uncertainty
    y_pred_all, std_all = predict_with_uncertainty(model, X_all_scaled.values)
    
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
    
    y_pred_all, _ = predict_with_uncertainty(model, X_all_scaled.values)
    
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
            'type': 'neural_network',
            'hidden_sizes': hidden_sizes,
            'dropout_rate': dropout_rate,
            'epochs': args.epochs,
            'learning_rate': learning_rate,
            'batch_size': batch_size,
            'patience': args.patience
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
        f.write(f"Cross-validation: {time.time() - cv_start:.2f} seconds\n")
        f.write(f"Predictions: {time.time() - pred_start:.2f} seconds\n")
        f.write(f"Metrics calculation: {time.time() - metrics_start:.2f} seconds\n")
        f.write(f"Feature importance calculation: {time.time() - importance_start:.2f} seconds\n")
        f.write(f"Results saving and plotting: {time.time() - save_start:.2f} seconds\n")

    # Save results to file if --save option is used
    if args.save:
        result_filename = f'{args.Y}_pred_nn_all_result.log'
        
        with open(result_filename, 'w') as f:
            f.write(f'Model: Neural Network\n')
            f.write(f'Features: {", ".join(args.X)}\n')
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

if __name__ == '__main__':
    print("Script started...")
    print("Importing required libraries...")
    main() 