import time
import optuna
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_validate, train_test_split
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic
from sklearn.ensemble import GradientBoostingRegressor as GBR
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.utils import shuffle

def main():
    parser = argparse.ArgumentParser(description='Perform linear regression using aggregated columns from multiple TSV files excluding the first column, calculate MAE, MSE, plot results, and save output.')
    parser.add_argument('--Y', required=True, help='File path for Y.tsv')
    parser.add_argument('--X', required=True, nargs='+', help='File paths for one or more X.tsv files')
    parser.add_argument('--C', default='concat_coord.tsv', help='File paths for one or more C.tsv files')
    parser.add_argument('--R', default='concat_row.tsv', help='File paths for one or more R.tsv files')
    parser.add_argument('--L', default='concat_element.tsv', help='File paths for one or more L.tsv files')
    parser.add_argument('-i', '--index', required=True, nargs='+', help='Column names to be used from the X.tsv files')
    parser.add_argument('-r', '--row', default=None, type=int)
    parser.add_argument('-c', '--coord', default=None, type=str)
    parser.add_argument('-z', '--zero', action='store_true', default=False)
    parser.add_argument('-o', '--output', dest='filename', type=str, default='', help="output filename")
    args = parser.parse_args()
    index = args.index
    row = args.row
    coord = args.coord
    zero = args.zero
    numb = len(index)
    
    filename = str(numb)
    if args.filename:
        filename = filename + '_' + args.filename
    if row:
        filename = filename + '_' + str(row) + 'd'
    if coord:
        filename = filename + '_' + coord
    if zero:
        filename = filename + '_zero'

    # Save results
    tsv_filename = f'gbr_optuna{filename}.tsv'
    png_filename = f'gbr_optuna{filename}.png'
    log_filename = f'gbr_optuna{filename}.log'
    
    # Load the data excluding the first column
    df_Y = pd.read_csv(args.Y, delimiter='\t').iloc[:, 1:]
    df_C = pd.read_csv(args.C, delimiter='\t', dtype=str).iloc[:, 1:]
    df_R = pd.read_csv(args.R, delimiter='\t', dtype=int).iloc[:, 1:]
    df_L = pd.read_csv(args.L, delimiter='\t', dtype=str).iloc[:, 1:]
    X_dataframes = []
    data_counts = []
    
    for x_file in args.X:
        df_X = pd.read_csv(x_file, delimiter='\t').iloc[:, 1:]
        melted_df = pd.melt(df_X)
        single_column_df = melted_df['value'].reset_index(drop=True)
        X_dataframes.append(single_column_df)
    
    df_X_combined = pd.concat(X_dataframes, axis=1)
    df_X_combined.columns = index
    df_Y_combined = pd.melt(df_Y.iloc[:df_X_combined.shape[0]])
    df_C_combined = pd.melt(df_C.iloc[:df_X_combined.shape[0]])
    df_R_combined = pd.melt(df_R.iloc[:df_X_combined.shape[0]])
    df_L_combined = pd.melt(df_L.iloc[:df_X_combined.shape[0]])
    
    X = df_X_combined
    Y = pd.DataFrame(df_Y_combined['value'])
    R = pd.DataFrame(df_R_combined['value'])
    L = pd.DataFrame(df_L_combined['value'])
    C = pd.DataFrame(df_C_combined['value'])
    
    Y.columns = ['Energy']
    R.columns = ['Row']
    L.columns = ['Metal']
    C.columns = ['Coordination']
    
    df_combined = pd.concat([R, L, C, X, Y], axis=1)
    df_combined = df_combined.dropna()
    # df_combined = df_combined[df_combined['Metal'] != 'Ba']
    if row:
        df_combined = df_combined[df_combined['Row'] == row]
    if coord:
        df_combined = df_combined[df_combined['Coordination'] == coord]
    if zero:
        df_combined = df_combined[df_combined['Energy'] < 0]
        
    X = df_combined.iloc[:, -(numb+1):-1]
    Y = df_combined['Energy']
    R = df_combined['Row']
    L = df_combined['Metal']
    C = df_combined['Coordination']

    # Split the data into training and test sets
    overall_start_time = time.time()
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

    # Define the objective function for Optuna
    def objective(trial):
        poly_degree = trial.suggest_int('poly__degree', 1, 3)
        n_estimators = trial.suggest_int('n_estimators', 50, 200)
        learning_rate = trial.suggest_float('learning_rate', 0.01, 0.1, log=True)
        subsample = trial.suggest_float('subsample', 0.8, 1.0)
        max_depth = trial.suggest_int('max_depth', 3, 5)
        min_samples_split = trial.suggest_int('min_samples_split', 2, 10)
        min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 4)
        max_features = trial.suggest_categorical('max_features', [None, 'sqrt', 'log2', 0.6, 0.8, 1.0])
        max_leaf_nodes = trial.suggest_categorical('max_leaf_nodes', [None, 10, 20, 30])
        min_weight_fraction_leaf = trial.suggest_float('min_weight_fraction_leaf', 0.0, 0.2)
        validation_fraction = trial.suggest_float('validation_fraction', 0.1, 0.2)
        n_iter_no_change = trial.suggest_int('n_iter_no_change', 10, 20)
        tol = trial.suggest_float('tol', 1e-4, 1e-2, log=True)
    
        # Create the pipeline with PolynomialFeatures, StandardScaler, and GradientBoostingRegressor
        gbr_pipe = Pipeline([
            ('poly', PolynomialFeatures(degree=poly_degree)),
            ('scaler', StandardScaler()),
            ('model', GBR(
                n_estimators=n_estimators,
                learning_rate=learning_rate,
                subsample=subsample,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                min_samples_leaf=min_samples_leaf,
                max_features=max_features,
                max_leaf_nodes=max_leaf_nodes,
                min_weight_fraction_leaf=min_weight_fraction_leaf,
                validation_fraction=validation_fraction,
                n_iter_no_change=n_iter_no_change,
                tol=tol,
                random_state=42
            )),
        ])
        
        # Cross-validate the pipeline and return the mean absolute error
        scores = cross_validate(gbr_pipe, X_train, Y_train, 
                                scoring='neg_mean_absolute_error', cv=5)
        mae = -np.mean(scores['test_score'])
        return mae

    # Create the Optuna study and optimize
    study = optuna.create_study(direction='minimize')
    start_time = time.time()
    study.optimize(objective, n_trials=1000)
    end_time = time.time()
    optimization_time = end_time - start_time
    
    # Extract the best parameters
    best_params = study.best_params
    best_poly_degree = best_params['poly__degree']
    best_n_estimators = best_params['n_estimators']
    best_learning_rate = best_params['learning_rate']
    best_subsample = best_params['subsample']
    best_max_depth = best_params['max_depth']
    best_min_samples_split = best_params['min_samples_split']
    best_min_samples_leaf = best_params['min_samples_leaf']
    best_max_features = best_params['max_features']
    best_max_leaf_nodes = best_params['max_leaf_nodes']
    best_min_weight_fraction_leaf = best_params['min_weight_fraction_leaf']
    best_validation_fraction = best_params['validation_fraction']
    best_n_iter_no_change = best_params['n_iter_no_change']
    best_tol = best_params['tol']

    # Create the best pipeline
    best_gbr_pipe = Pipeline([
        ('poly', PolynomialFeatures(degree=best_poly_degree)),
        ('scaler', StandardScaler()),
        ('pca', PCA(n_components=0.95)),  # Keep 95% variance
        ('model', GBR(
            n_estimators=best_n_estimators,
            learning_rate=best_learning_rate,
            subsample=best_subsample,
            max_depth=best_max_depth,
            min_samples_split=best_min_samples_split,
            min_samples_leaf=best_min_samples_leaf,
            max_features=best_max_features,
            max_leaf_nodes=best_max_leaf_nodes,
            min_weight_fraction_leaf=best_min_weight_fraction_leaf,
            validation_fraction=best_validation_fraction,
            n_iter_no_change=best_n_iter_no_change,
            tol=best_tol,
            random_state=42
        )),
    ])
    
    # Fit the best pipeline to the training data
    start_time = time.time()
    best_gbr_pipe.fit(X_train, Y_train)
    end_time = time.time()
    fitting_time = end_time - start_time
    
    # Log the best parameters and consumed time
    with open(log_filename, 'w') as file:
        file.write(f"Optimized poly: {best_poly_degree}\n")
        file.write(f"Optimized n_estimators: {best_n_estimators}\n")
        file.write(f"Optimized learning_rate: {best_learning_rate:.4f}\n")
        file.write(f"Optimized subsample: {best_subsample:.4f}\n")
        file.write(f"Optimized max_depth: {best_max_depth}\n")
        file.write(f"Optimized min_samples_split: {best_min_samples_split}\n")
        file.write(f"Optimized min_samples_leaf: {best_min_samples_leaf}\n")
        file.write(f"Optimized max_features: {best_max_features}\n")
        file.write(f"Optimized max_leaf_nodes: {best_max_leaf_nodes}\n")
        file.write(f"Optimized min_weight_fraction_leaf: {best_min_weight_fraction_leaf:.4f}\n")
        file.write(f"Optimized validation_fraction: {best_validation_fraction:.4f}\n")
        file.write(f"Optimized n_iter_no_change: {best_n_iter_no_change}\n")
        file.write(f"Optimized tol: {best_tol:.4f}\n\n")

    # Cross-validate the pipeline and print CV scores for GBR
    start_time = time.time()
    gbr_score = cross_validate(best_gbr_pipe, X_train, Y_train, 
                               scoring=['r2', 'neg_mean_absolute_error', 'neg_mean_squared_error'], cv=5)
    end_time = time.time()
    cross_validation_time = end_time - start_time
    
    with open(log_filename, 'a') as file:
        file.write(f"\tR^2\tMAE\tMSE\n")
        file.write(f"CV\t{np.mean(gbr_score['test_r2']):.4f}\t{-np.mean(gbr_score['test_neg_mean_absolute_error']):.4f}\t{-np.mean(gbr_score['test_neg_mean_squared_error']):.4f}\n")
    
    # Predict on the entire set using the final GBR model
    start_time = time.time()
    Y_pred_gbr = best_gbr_pipe.predict(X)
    end_time = time.time()
    prediction_time = end_time - start_time
    
    # Compute and print MAE and MSE for the entire set for GBR
    mae_gbr = mean_absolute_error(Y, Y_pred_gbr)
    mse_gbr = mean_squared_error(Y, Y_pred_gbr)
    
    with open(log_filename, 'a') as file:
        file.write(f"Entire\t{best_gbr_pipe.score(X, Y):.4f}\t{mae_gbr:.4f}\t{mse_gbr:.4f}\n")
    
    # Predict on the test set using the final GBR model
    start_time = time.time()
    Y_pred_gbr_test = best_gbr_pipe.predict(X_test)
    end_time = time.time()
    test_prediction_time = end_time - start_time
    
    # Compute and print MAE and MSE for the test set for GBR
    mae_gbr_test = mean_absolute_error(Y_test, Y_pred_gbr_test)
    mse_gbr_test = mean_squared_error(Y_test, Y_pred_gbr_test)
    
    with open(log_filename, 'a') as file:
        file.write(f"Test\t{best_gbr_pipe.score(X_test, Y_test):.4f}\t{mae_gbr_test:.4f}\t{mse_gbr_test:.4f}\n\n")
    
    overall_end_time = time.time()
    overall_time = overall_end_time - overall_start_time
    
    with open(log_filename, 'a') as file:
        file.write(f"Optimization time: {optimization_time:.2f} sec\n")
        file.write(f"Model fitting time: {fitting_time:.2f} sec\n")
        file.write(f"Cross-validation time: {cross_validation_time:.2f} sec\n")
        file.write(f"Prediction time (entire set): {prediction_time:.2f} sec\n")
        file.write(f"Prediction time (test set): {test_prediction_time:.2f} sec\n")
        file.write(f"Overall time: {overall_time:.2f} sec\n")
    
    # Save predictions and residuals to a TSV file
    df_combined['Predicted Energy'] = Y_pred_gbr
    df_combined['Residuals'] = Y - Y_pred_gbr
    df_combined.to_csv(tsv_filename, sep='\t', index=False)

#     # Plot results
#     plt.figure(figsize=(10, 8))
#     colors = ['red', 'green', 'blue']
#     markers = ['v', '^', 's', 'D', 'o']
#     for i, row in enumerate([3, 4, 5]):
#         sub = df_combined[df_combined['Row'] == row]
#         for j, coordination in enumerate(['WZ', 'ZB', 'TN', '33', 'RS']):
#             subset = sub[sub['Coordination'] == coordination]
#             plt.scatter(subset['Energy'], subset['Predicted Energy'], alpha=0.3, color=colors[i], marker=markers[j], label=f'{row}_{coordination}')
#             for (x, y, label) in zip(subset['Energy'], subset['Predicted Energy'], subset['Metal']):
#                 plt.annotate(label, (x, y))
    
#     plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], 'r--', lw=1)
#     plt.xlabel('DFT-calculated Formation Energy (eV)')
#     plt.ylabel('Predicted Formation Energy (eV)')
#     plt.legend()
#     plt.tight_layout()
#     plt.gcf().savefig(png_filename, bbox_inches="tight")
#     plt.close()

#     # Plot correlation matrix
#     M = pd.concat([Y, X], axis=1)
#     correlation_matrix = M.corr()
#     abs_correlation_matrix = correlation_matrix.abs()
    
#     plt.figure(figsize=(7, 6))
#     sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='coolwarm')
#     plt.xticks(np.arange(M.shape[1]) + 0.5, M.columns, rotation=90, ha='right')
#     plt.yticks(np.arange(M.shape[1]) + 0.5, M.columns, rotation=0, va='center')
#     plt.tight_layout()
#     plt.gcf().savefig(f'covariance_matrix{filename}.png', bbox_inches="tight")
#     plt.close()
    
#     plt.figure(figsize=(7, 6))
#     sns.heatmap(abs_correlation_matrix, annot=True, fmt=".2f", cmap='coolwarm', center=0, vmin=0, vmax=1)
#     plt.xticks(np.arange(M.shape[1]) + 0.5, M.columns, rotation=90, ha='right')
#     plt.yticks(np.arange(M.shape[1]) + 0.5, M.columns, rotation=0, va='center')
#     plt.tight_layout()
#     plt.gcf().savefig(f'abs_covariance_matrix{filename}.png', bbox_inches="tight")
#     plt.close()

if __name__ == "__main__":
    main()
