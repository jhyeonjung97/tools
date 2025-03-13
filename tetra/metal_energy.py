import os
import numpy as np
import pandas as pd

root = '/pscratch/sd/j/jiuy97'

kjmol = 96.485
calmol = 23.061

h2o = -14.23091949
cvh2o = 0.103
tsh2o = 0.675
zpeh2o = 0.558
gh2o = h2o + cvh2o - tsh2o + zpeh2o

h2 = -6.77149190
cvh2 = 0.0905
tsh2 = 0.408
zpeh2 = 0.268
gh2 = h2 + cvh2 - tsh2 + zpeh2
gh = gh2 / 2

go = gh2o - gh2 + 2.46
ho = -5.9248399 / 2

metals = [
    # metal, #M, #O, #N, G_form, E_oxide, Gcorr_oxide
    ['Ti', 1, 2, 4, -212.330, -99.45185983, 0.548526],
    ['V',  2, 5, 2, -344.000, -104.1935893, 0.673172],
    ['Cr', 2, 3, 2, -260.530, -79.47591331, 0.621017],
    ['Mn', 1, 1, 2, -90.210,  -31.17185317, 0.046701],
    ['Fe', 2, 3, 2, -177.100, -68.25294948, 0.449591],
    ['Co', 3, 4, 2, -167.835, -83.73529149, 0.572332],
    ['Ni', 1, 1, 2, -51.610,  -20.40998666, 0.124372],
    ['Cu', 1, 1, 2, -30.400,  -17.83506218, 0.140476],
]

columns = ['Metal', '#M', '#O', '#N', 'G_form', 'E_oxide', 'Gcorr_oxide']
df = pd.DataFrame(metals, columns=columns)

df['G_form'] = df['G_form'] / df['#M'] / calmol
df['G_oxide'] = (df['E_oxide'] + df['Gcorr_oxide']) / df['#M'] / df['#N']
df['G_metal'] = df['G_oxide'] - df['G_form'] - df['#O'] / df['#M'] * go

print(df)