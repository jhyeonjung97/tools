import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Perform linear regression using aggregated columns from multiple TSV files excluding the first column, calculate MAE, MSE, plot results, and save output.')
    parser.add_argument('--X', required=True, nargs='+', help='File paths for one or more X.tsv files')
    parser.add_argument('-o', '--output', dest='filename', type=str, default='', help="output filename")
    args = parser.parse_args()
    filename = args.filename

    X_dataframes = []
    for x_file in args.X:
        df_X = pd.read_csv(x_file, delimiter='\t').iloc[:, 1:]
        X_dataframes.append(df_X)
        
    df_X_combined = pd.concat(X_dataframes, axis=0)
    # df_X_combined = df_X_combined.astype(float)
    df_X_combined.to_csv(f'concat_{filename}.tsv', sep='\t', index=True, float_format='%.2f')

if __name__ == "__main__":
    main()
