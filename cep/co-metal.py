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
    # Name/Ueff, #M, #O, #H, #N, G_form, E_oxide, Gcorr_oxide
    ['CoLH', 1, 2, 2, 4, -109.999, -107.83367123, 0.0],
]

columns = ['Metal', '#M', '#O', '#H', '#N', 'G_form', 'E_oxide', 'Gcorr_oxide']
df = pd.DataFrame(metals, columns=columns)

df['G_form'] = df['G_form'] / df['#M'] / calmol
df['G_oxide'] = (df['E_oxide'] + df['Gcorr_oxide'] ) / df['#M'] / df['#N']
df['G_metal'] = df['G_oxide'] - df['G_form'] - df['#O'] / df['#M'] * go - df['#H'] / df['#M'] * gh

print(df)