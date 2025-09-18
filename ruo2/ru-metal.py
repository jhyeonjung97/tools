import os
import numpy as np
import pandas as pd

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
    # Ueff, #M, #O, #N, G_form, E_oxide, Gcorr_oxide
    [0, 1, 2, 8, -40.700, -178.20963279, 0.0],
    [1, 1, 2, 8, -40.700, -170.09190453, 0.0],
    [2, 1, 2, 8, -40.700, -162.78761804, 0.0],
    [3, 1, 2, 8, -40.700, -156.32933730, 0.0],
    [4, 1, 2, 8, -40.700, -150.34186082, 0.0],
]

Gcorr_oxide = -9.243204045 + 2*go -40.700/calmol +178.20963279 / 8

columns = ['Metal', '#M', '#O', '#N', 'G_form', 'E_oxide', 'Gcorr_oxide']
df = pd.DataFrame(metals, columns=columns)

df['G_form'] = df['G_form'] / df['#M'] / calmol
df['G_oxide'] = (df['E_oxide'] + df['Gcorr_oxide'] ) / df['#M'] / df['#N'] #+ Gcorr_oxide
df['G_metal'] = df['G_oxide'] - df['G_form'] - df['#O'] / df['#M'] * go

print(df)