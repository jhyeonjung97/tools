# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2

calmol = 23.061

mno2 = -190.724224/8
# mn = -13.303998/2 # -6.651999
mn = -.52017828E+03/58
formation_energy = -111.100/calmol

mn_ueff = mno2 - formation_energy - 2 * gh2o + 2 * gh2

print(mn,mn_ueff)