import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from matplotlib.font_manager import FontProperties

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
    
    # Load the data excluding the first column
    df_Y = pd.read_csv(args.Y, delimiter='\t').iloc[:, 1:]
    df_C = pd.read_csv(args.C, delimiter='\t', dtype=str).iloc[:, 1:]
    df_R = pd.read_csv(args.R, delimiter='\t', dtype=int).iloc[:, 1:]
    df_L = pd.read_csv(args.L, delimiter='\t', dtype=str).iloc[:, 1:]
    X_dataframes = []
    
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
    
    Y.columns = ['Formation Energy']
    R.columns = ['Row']
    L.columns = ['Metal']
    C.columns = ['Coordination']
    
    df_combined = pd.concat([R, L, C, X, Y], axis=1)
    df_combined = df_combined.dropna()

    if row:
        df_combined = df_combined[df_combined['Row'] == row]
    if coord:
        df_combined = df_combined[df_combined['Coordination'] == coord]
    if zero:
        df_combined = df_combined[df_combined['Formation Energy'] < 0]
        
    X = df_combined.iloc[:, -(numb+1):-1]
    Y = df_combined['Formation Energy']
    R = df_combined['Row']
    L = df_combined['Metal']
    C = df_combined['Coordination']

    model = LinearRegression()
    model.fit(X, Y)

    Y_pred = model.predict(X)
    Y_pred = pd.Series(Y_pred, index=Y.index)  # Convert Y_pred to a pandas Series to match Y
    
    mae = mean_absolute_error(Y, Y_pred)
    mse = mean_squared_error(Y, Y_pred)

    tsv_filename = f'regression{filename}.tsv'
    png_filename = f'regression{filename}.png'
    log_filename = f'regression{filename}.log'
    
    with open(log_filename, 'w') as file:
        file.write(f"\nIntercept: {model.intercept_}\n\n")
        for i, coef in enumerate(model.coef_):
            file.write(f"Coefficient ({index[i]}): {coef:.2f}\n")
        file.write(f"\nR-squared: {model.score(X, Y):.4f}\n")
        file.write(f"Mean Absolute Error: {mae:.4f}\n")
        file.write(f"Mean Squared Error: {mse:.4f}\n\n")
    
    df_combined['Predicted'] = Y_pred
    df_combined['Residuals'] = Y - Y_pred
    df_combined.to_csv(tsv_filename, sep='\t', index=False)
        
    plt.figure(figsize=(8, 6), dpi=300)
    colors = ['red', 'green', 'blue']
    markers = ['>', '<', 'o', 's', 'p', 'd', 'h', '^', 'v']
    for i, row in enumerate([3, 4, 5]):
        sub = df_combined[df_combined['Row'] == row]
        for j, coordination in enumerate(['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT', 'AQ', 'AU']):
            subset = sub[sub['Coordination'] == coordination]
            LL = subset['Metal']
            YY = subset['Formation Energy']
            YY_pred = subset['Predicted']
            plt.scatter(YY, YY_pred, alpha=0.3, color=colors[i], marker=markers[j], label=f'{row}_{coordination}')
            for (x, y, label) in zip(YY, YY_pred, LL):
                plt.annotate(label, (x, y), fontsize=8)
            
    plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], 'r--', lw=1)
    plt.xlabel('DFT-calculated Energy (eV)')
    plt.ylabel('Predicted Energy (eV)')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), handlelength=1.0, labelspacing=0.3)
    plt.tight_layout()
    plt.savefig(png_filename, bbox_inches="tight")
    plt.close()
    print(f"Figure saved as {png_filename}")

    M = pd.concat([Y, X], axis=1)
    covariance_matrix = M.cov()
    correlation_matrix = M.corr()
    
    covariance_matrix_filename = f'covariance_matrix{str(filename)}.tsv'
    covariance_matrix.to_csv(covariance_matrix_filename, sep='\t')
    
    plt.figure(figsize=(8, 6), dpi=300) # Set the figure size as needed
    ax = sns.heatmap(covariance_matrix, annot=True, fmt=".2f", annot_kws={"size": 5}, cmap='coolwarm', cbar=False)
    ax.set_xticks(np.arange(M.shape[1]) + 0.5)
    ax.set_xticklabels(M.columns, rotation=90, ha='right', fontsize=8)
    ax.set_yticks(np.arange(M.shape[1]) + 0.5)
    ax.set_yticklabels(M.columns, rotation=0, va='center', fontsize=8)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    cbar = plt.colorbar(ax.collections[0], ax=ax, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_visible(False)
    plt.tight_layout()
    plt.savefig(f'covariance_matrix{str(filename)}.png', bbox_inches="tight")
    plt.close()

    correlation_matrix_filename = f'correlation_matrix{str(filename)}.tsv'
    correlation_matrix.to_csv(correlation_matrix_filename, sep='\t')
    
    plt.figure(figsize=(8, 6), dpi=300) # Set the figure size as needed
    ax = sns.heatmap(correlation_matrix, annot=True, fmt=".2f", annot_kws={"size": 5}, cmap='coolwarm', cbar=False)
    # colors = ['black', 
    #           '#FFCCCC',  # light-red
    #           '#FFCCCC',  # light-red
    #           'red', 
    #           'black',
    #           '#FFCCCC',  # light-red
    #           '#FFCCCC',  # light-red
    #           '#81C784',  # darker light-yellow
    #           '#81C784',  # darker light-yellow
    #           'black',
    #           'grey',
    #           '#FFEB3B',  # softer light-green
    #           '#81C784',  # darker light-yellow
    #           '#FFEB3B',  # softer light-green
    #           '#81C784',  # darker light-yellow
    #           '#81C784',  # darker light-yellow
    #           '#81C784',  # darker light-yellow
    #           '#81C784',  # darker light-yellow
    #           '#FFCCCC',  # light-red
    #           '#81C784',  # darker light-yellow
    #           'green',  # yellow-orange
    #           '#FFC107',
    #           '#81C784',  # darker light-yellow
    #           '#81C784',  # darker light-yellow
    #           'green',  # yellow-orange
    #           '#64B5F6',  # softer light-blue
    #           '#64B5F6',  # softer light-blue
    #           '#64B5F6',  # softer light-blue
    #           'blue',
    #           '#64B5F6',  # softer light-blue
    #           '#64B5F6',  # softer light-blue
    #          ]
    # light_colors_indices = [1, 2, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 25, 26, 27, 29, 30]
    # bold_indices = [i for i, color in enumerate(colors) if color != 'grey' and i not in light_colors_indices]
    ax.set_xticks(np.arange(M.shape[1]) + 0.5)
    ax.set_xticklabels(M.columns, rotation=90, ha='right', fontsize=8)
    ax.set_yticks(np.arange(M.shape[1]) + 0.5)
    # for i, (tick_label, color) in enumerate(zip(ax.get_yticklabels(), colors)):
    #     if i in bold_indices:
    #         tick_label.set_fontproperties(FontProperties(weight='bold'))
    #     tick_label.set_color(color)
    ax.set_yticklabels(M.columns, rotation=0, va='center', fontsize=8)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    cbar = plt.colorbar(ax.collections[0], ax=ax, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_visible(False)
    plt.tight_layout()
    plt.savefig(f'correlation_matrix{str(filename)}.png', bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
