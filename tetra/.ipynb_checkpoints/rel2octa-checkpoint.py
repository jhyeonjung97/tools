import pandas as pd
import numpy as np
import sys

def read(data):
    return pd.read_csv(data, delimiter='\t')

def write(df, filename):
    df.to_csv(filename, sep='\t', index=False)

data = sys.argv[1]
df = read(data).iloc[:, 1:]
filename, ext = data.rsplit('.', 1)

for metal_row in ['3d', '4d', '5d']:
    for i in range(6):
        for j in range(13):
            index1 = 13 * i + j
            index2 = 13 * 5 + j
            if index1 < len(df) and index2 < len(df):
                if not pd.isna(df.at[index1, metal_row]) and not pd.isna(df.at[index2, metal_row]):
                    df.at[index1, metal_row] = df.at[index1, metal_row] - df.at[index2, metal_row]

df = df[:-13]
index_pattern = np.tile(np.arange(13), len(df) // 13 + 1)[:len(df)]
df.index = index_pattern
df.to_csv(f'{filename}_rel.{ext}', sep='\t', index=True)
print(f"DataFrame saved as {filename}_rel.{ext}")