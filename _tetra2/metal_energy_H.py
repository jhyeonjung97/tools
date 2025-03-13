import os
import numpy as np
import pandas as pd

root = '/pscratch/sd/j/jiuy97'

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

metal_rows = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']
}

columns = ['Metal', '#M', '#O', 'G_form', 'N_oxide', 'E_oxide', 'Gcorr_oxide', 'N_metal', 'E_metal']

metals = [
    # metal, #M, #O, G_form, N_oxide, E_oxide, Gcorr_oxide, N_metal, E_metal
    ['Ti', 1, 2, -212.330, 4, -99.45185983, 0.548526, 2, -11.2546574],
    ['V',  2, 5, -344.000, 2, -104.1935893, 0.673172, 2, -11.9316278],
    ['Cr', 2, 3, -260.530, 2, -79.47591331, 0.621017, 2, -13.8032561],
    ['Mn', 1, 1, -90.210,  2, -31.17185317, 0.046701, 2, -13.3039980],
    ['Fe', 2, 3, -177.100, 2, -68.25294948, 0.449591, 2, -9.58156656],
    ['Co', 3, 4, -167.835, 2, -83.73529149, 0.572332, 2, -8.71118563],
    ['Ni', 1, 1, -51.610,  2, -20.40998666, 0.124372, 2, -3.14872124],
    ['Cu', 1, 1, -30.400,  2, -17.83506218, 0.140476, 2, -5.33647942],
]

df = pd.DataFrame(metals, columns=columns)

df['G_oxide'] = (df['E_oxide'] + df['Gcorr_oxide']) / df['N_oxide']
# df['G_metal'] = (df['E_metal'] + df['Gcorr_metal']) / df['N_metal']
# df['G_metal'] = df['E_metal'] / df['N_metal']
df['G_metal'] = df['G_oxide'] - df['G_form'] - df['#O'] / df['#M'] * go

print(df)