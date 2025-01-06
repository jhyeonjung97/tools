# Water and hydrogen properties
water_E = -14.23797429
water_Cv = 0.103
water_TS = 0.675
water_ZPE = 0.558
water_G = water_E + water_Cv - water_TS + water_ZPE

hydrogen2_E = -6.77412273
hydrogen2_Cv = 0.0905
hydrogen2_TS = 0.408
hydrogen2_ZPE = 0.273
hydrogen2_G = hydrogen2_E + hydrogen2_Cv - hydrogen2_TS + hydrogen2_ZPE
hydrogen_G = hydrogen2_G / 2

oxygen_G = water_G - hydrogen2_G
hydroxide_G = water_G - hydrogen_G

# Correction values
OH_Cv = 0.042
OH_TS = 0.066
OH_ZPE = 0.376
OH_corr = OH_Cv - OH_TS + OH_ZPE

O_Cv = 0.034
O_TS = 0.060
O_ZPE = 0.064
O_corr = O_Cv - O_TS + O_ZPE

OOH_Cv = 0.077
OOH_TS = 0.134
OOH_ZPE = 0.471
OOH_corr = OOH_Cv - OOH_TS + OOH_ZPE

E_ = -279.1
E_OH = -289.57160323
E_O = -284.14645182
E_OOH = -294.13407435

G_ = E_
G_OH = E_OH + OH_corr
G_O = E_O + O_corr
G_OOH = E_OOH + OOH_corr

dG_OH = G_OH - G_ - hydroxide_G
dG_O = G_O - G_ - oxygen_G
dG_OOH = G_OOH - G_ - oxygen_G - hydroxide_G

dG1 = dG_OH
dG2 = dG_O - dG_OH
dG3 = dG_OOH - dG_O
dG4 = 4.92 - dG_OOH

OER = max(dG1, dG2, dG3, dG4) - 1.23
ORR = 1.23 - min(dG1, dG2, dG3, dG4)

print(G_, G_OH, G_O, G_OOH, dG_OH, dG_O, dG_OOH, dG1, dG2, dG3, dG4, OER, ORR)





















