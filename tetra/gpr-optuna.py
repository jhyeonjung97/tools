import time
import optuna
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV, cross_validate, train_test_split
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.linear_model import Ridge
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
    tsv_filename = f'gpr_optuna{filename}.tsv'
    png_filename = f'gpr_optuna{filename}.png'
    log_filename = f'gpr_optuna{filename}.log'
            
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
        alpha = trial.suggest_float('model__alpha', 1e-10, 10, log=True)
    
        # Create the pipeline with PolynomialFeatures and StandardScaler for GPR
        gpr_pipe = Pipeline([
            ('poly', PolynomialFeatures(degree=poly_degree)),
            ('scaler', StandardScaler()),
            ('model', GPR(alpha=alpha, normalize_y=True)),
        ])
    
        # Cross-validate the pipeline and return the mean absolute error
        gpr_score = cross_validate(gpr_pipe, X_train, Y_train, 
                                   scoring='neg_mean_absolute_error', cv=5)
        mae = -np.mean(gpr_score['test_score'])
        return mae
    
    # Create the Optuna study and optimize
    start_time = time.time()
    study = optuna.create_study(direction='minimize')
    study.optimize(objective, n_trials=100)
    end_time = time.time()
    optimization_time = end_time - start_time

    # Extract the best parameters and pipeline
    best_params = study.best_params
    best_poly_degree = best_params['poly__degree']
    best_alpha = best_params['model__alpha']
    
    # Create the best pipeline
    best_gpr_pipe = Pipeline([
        ('poly', PolynomialFeatures(degree=best_poly_degree)),
        ('scaler', StandardScaler()),
        ('model', GPR(alpha=best_alpha, normalize_y=True)),
    ])
    
    # Fit the best pipeline to the training data
    start_time = time.time()
    best_gpr_pipe.fit(X_train, Y_train)
    end_time = time.time()
    fitting_time = end_time - start_time
    
    # Log the best parameters
    with open(log_filename, 'w') as file:
        file.write(f"Optimized poly: {best_poly_degree}\n")
        file.write(f"Optimized alpha: {best_alpha:.4f}\n\n")
    
    # Cross-validate the pipeline and print CV scores for GPR
    start_time = time.time()
    gpr_score = cross_validate(best_gpr_pipe, X_train, Y_train, 
                               scoring=['r2', 'neg_mean_absolute_error', 'neg_mean_squared_error'], cv=5)
    end_time = time.time()
    cross_validation_time = end_time - start_time
    
    with open(log_filename, 'a') as file:
        file.write("\tR^2\tMAE\tMSE\n")
        file.write(f"CV\t{np.mean(gpr_score['test_r2']):.4f}\t{-np.mean(gpr_score['test_neg_mean_absolute_error']):.4f}\t{-np.mean(gpr_score['test_neg_mean_squared_error']):.4f}\n")
    
    # Predict on the entire set using the final GPR model
    start_time = time.time()
    Y_pred_gpr = best_gpr_pipe.predict(X)
    end_time = time.time()
    prediction_time = end_time - start_time
    
    # Compute and print MAE and MSE for the entire set for GPR
    mae_gpr = mean_absolute_error(Y, Y_pred_gpr)
    mse_gpr = mean_squared_error(Y, Y_pred_gpr)
    
    with open(log_filename, 'a') as file:
        file.write(f"Entire\t{best_gpr_pipe.score(X, Y):.4f}\t{mae_gpr:.4f}\t{mse_gpr:.4f}\n")
    
    # Predict on the test set using the final GPR model
    start_time = time.time()
    Y_pred_gpr_test = best_gpr_pipe.predict(X_test)
    end_time = time.time()
    test_prediction_time = end_time - start_time
    
    # Compute and print MAE and MSE for the test set for GPR
    mae_gpr_test = mean_absolute_error(Y_test, Y_pred_gpr_test)
    mse_gpr_test = mean_squared_error(Y_test, Y_pred_gpr_test)
    
    with open(log_filename, 'a') as file:
        file.write(f"Test\t{best_gpr_pipe.score(X_test, Y_test):.4f}\t{mae_gpr_test:.4f}\t{mse_gpr_test:.4f}\n\n")
        
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
    df_combined['Predicted Energy'] = Y_pred_gpr
    df_combined['Residuals'] = Y - Y_pred_gpr
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
