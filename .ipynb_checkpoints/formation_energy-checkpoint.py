import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

file_path = '/pscratch/sd/j/jiuy97/3_V_shape/metal/merged_norm_energy.tsv'  # Adjust the file extension if necessary
df = pd.read_csv(file_path, delimiter='\t').iloc[:, 1:]

# Select the first three columns
first_three_columns = df.iloc[:, :3]

# Compute the minimum of these columns for each row
min_values = first_three_columns.min(axis=1)

# If you want to replace these three columns with the minimum values column
df = df.drop(df.columns[:3], axis=1)  # This removes the first three columns
df.insert(0, 'Min_of_first_three', min_values)  # This inserts the min_values column as the first column

print(df)
# Or, if you want to create a new DataFrame
# new_df = pd.DataFrame(min_values, columns=['Min_of_first_three'])
# new_df = pd.concat([new_df, df.iloc[:, 3:]], axis=1)  # Add the rest of the columns from the original df

# new_df now contains the minimum values of the first three columns and the rest of the original columns