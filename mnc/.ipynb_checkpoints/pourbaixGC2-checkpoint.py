### Evan Carlson ###
### Date: 2024-02-25 ###
### This script generates surface pourbaix diagrams from Grand Canonical DFT calculations ###

### Import Dependencies ###
from math import pow
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.ticker import *
from pylab import *
from scipy import *

### Define New Figure ###
figName = 'SurfacePBX_Test'
fig = plt.figure(figsize=(6.5,6),dpi=300)
ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])
  
# Voltage Range for PBX Diagram
Emin = -0.9
Emax = 2.1

# Values for X,Y Limits of the Plot
pH = np.arange(0, 14.1, 0.01)
E = np.arange(Emin, Emax, 0.01) 

# Set Axes Limits and Labels
ax.axis([0, 14, Emin, Emax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

### Constants ###
kbt = 0.0256  # RT/nF @ 25C
RHE_Shift = kbt*np.log(10)  # 59mV
kjmol = 96.485  # kJ/eV


##### Possible Coverages in Surface PBX #####
## First entry has to be pristine surface (all other entries are relative to this one)
## Entries should be added in the following format: 
## [Vacuum DFT Energy, #Hs, #Os, #OHs, #OOHs, #H2Os, A, B, C]
## #Hs, #Os, #OHs, #OOHs, #H2Os, are the number of H, O, OH, OOH, H2O adsorbates on the surface relative to the first entry
## A, B, C are the coefficients of a parabolic fit (G = A*U^2 + B*U + C) of the Grand Canonical DFT Energy (U is the applied potential vs. SHE)

surfs = [    
    [-420.51205470, 2, 0, 0, 0, -0.33415314945322344, -0.10794642676914794, -421.8598852435611],  # OBridge_vacTop (Rutile 110)
    [-439.61985055, 0, 0, 0, 0, -0.36597110528046517, 0.0664778573230075, -441.19851466871967],  # OBridge_OHTop_fromWater (Rutile 110)
    [-427.68085946, 1, 0, 0, 0, -0.3714150738091594, 0.21734914901391614, -429.63978276116586],  # OBridge_OTop (Rutile 110)
    [-451.33981839, 1, 0, 0, 0, -0.3438412618178991, -0.33615544648583073, -452.2326262543747],  # OHBridge_OHTopSide2 (Rutile 110)
    [-431.10975082, 2, 0, 0, 0, -0.36083879711605205, -0.6948350818248356, -431.93789281079404],  # OHBridge_vacTop2 (Rutile 110)
    [-588.42781682, 0, 0, 0, 0, -0.8299846678733894, -0.4206383098177547, -592.8773643194726],  # OBridge_vacTop (Alpha 110)
    [-627.79024904, 0, 0, 4, 0, 0, -0.7036194647471677, 0.3161797244963896, -633.5198738630697],  # OBridge_OHTop (Alpha 110)
    [-604.61945692, 0, 4, 0, 0, 0, -0.7214245472257269, 0.31518120220547347, -610.2183732425576],  # OBridge_OTop (Alpha 110)
    [-604.06827920, 3, 0, 0, 0, 0, -0.9680697130706117, -0.9980461956425379, -607.612927049476],  # OHBridgeABC_OBridgeD_VacTop (Alpha 110)
    [-608.72837309, 4, 0, 0, 0, 0, -0.9970808111805852, -1.1224597833278684, -612.2014148705357],  # OHBridge_vacTop (Alpha 110)
    [-627.59316584, 4, 4, 0, 0, 0, -0.8133854726984423, 0.060545828520655695, -632.3953466194389],  # OHBridge_OTop (Alpha 110)
    [-649.96329741, 4, 0, 4, 0, 0, -0.9296934058838175, -0.7490523053555069, -653.2410999197407],  # OHBridge_OHTop (Alpha 110)
    [-588.42781682, 0, 0, 0, 0, 0, -0.8299846678733894, -0.4206383098177547, -592.8773643194726],  # Dummy OBridge_vacTop (Alpha 110)
    [-645.21429491, 3, 0, 4, 0, 0, -0.838476241591988, -0.6094077071550351, -648.7112369556921],  # OHBridgeABC_OBridgeD_OHTop (Alpha 110)
    [-621.98399886, 3, 4, 0, 0, 0, -0.7099814621827858, 0.16383757179439587, -627.2593072067509],  # OHBridgeABC_OBridgeD_OTop (Alpha 110)
    [-438.12892127, 0, 0, 0, 0, 0, -0.7246730196158068, 1.7570034771953849, -441.1449285616018],  # OBridge_vacTop (Alpha 100)
    [-458.14880201, 0, 0, 2, 0, 0, -0.6535615161502478, 2.109495758112015, -461.8944013962121],  # OBridge_OHTop (Alpha 100)
    [-446.45300783, 0, 2, 0, 0, 0, -0.6428986913146885, 2.1679290355605363, -450.44212231004127],  # OBridge_OTop (Alpha 100)
    [-468.9081545, 2, 0, 2, 0, 0, -0.691786023339257, 1.3096963819989538, -470.90773609497387],  # OHBridge_OHTop (Alpha 100)
    [-448.13516333, 2, 0, 0, 0, 0, -0.7158333306459814, 0.8580121098672541, -450.2448435169503],  # OHBridge_vacTop (Alpha 100)
    [-458.05287266, 1, 1, 1, 0, 0, -0.7384235923910929, 2.0031699089410298, -461.19596535353253],  # OHBridgeA_OBridgeB_O+OHTop (relaxed OHBridge_OTop) (Alpha 100)
    [-443.68338301, 1, 0, 0, 0, 0, -0.6787628442455168, 1.2122923154925143, -445.90509038918333],  # OHBridgeA_OBridgeB_VacTop **VRANGE_POS** (Alpha 100)
    [-443.68338301, 1, 0, 0, 0, 0, -0.625789707123459, 1.059679892435809, -445.82804530248416],  # OHBridgeA_OBridgeB_VacTop **FULL**(Alpha 100)
    [-443.68338301, 1, 0, 0, 0, 0, -0.563223073798204, 1.0037196190764035, -445.8658128549143],  # OHBridgeA_OBridgeB_VacTop **VRANGE_NEG**(Alpha 100)
    [-438.12892127, 0, 0, 0, 0, 0, -0.7246730196158068, 1.7570034771953849, -441.1449285616018],  # Dummy OBridge_vacTop (Alpha 100)    
    [-464.05935517, 1, 0, 2, 0, 0, -0.5941122121140239, 1.3796041244437771, -466.40497043739225],  # OHBridgeA_OBridgeB_OHTop (Alpha 100)
    [-452.30535743, 1, 2, 0, 0, 0, -0.7731421167692354, 2.2636758482510673, -455.8063515131638],  # OHBridgeA_OBridgeB_OTop (Alpha 100)
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-446.75018392, 0, 0, 2, 0, 0, -0.86407420024696, 3.4779271259934883, -452.8218029284715],  # OBridge_OHTop (Alpha (100))
    [-434.84285041, 0, 2, 0, 0, 0, -0.8345538999405395, 3.493379228883208, -441.16221871033605],  # OBridge_OTop (Alpha (100))
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-437.6968926, 2, 0, 0, 0, 0, -0.6684369486354754, 1.39726870907811, -441.06365540327977],  # OHBridge_VacTop (Alpha (100))
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-458.53866205, 2, 0, 2, 0, 0, -0.5982872194885003, 1.6208817186937816, -461.65564909392157],  # OHBridge_OHTop (Alpha (100) ****VRange****        
    [-427.15021798, 0, 0, 0, 0, 0, -0.8513001772902745, 2.5589486132301813, -431.6800432612759],  # OBridge_vacTop (Alpha (100))
    [-578.6959023, 0, 0, 0, 0, 0, -1.2769268721820428, 0.37165879352215103, -583.9474528155971],  # OBridge_vacTop (Alpha (110))
    [-617.30877079, 0, 0, 4, 0, 0, -0.7385894503552835, 1.23156625746479, -624.6491911419675],  # OBridge_OHTop (Alpha (110))
    [-593.70496404, 0, 4, 0, 0, 0, -0.6579588283505265, 1.0231986986924992, -601.0027785401438],  # OBridge_OTop (Alpha (110)
    [-640.49878946, 4, 4, 0, 0, 0, -0.8469360154918735, -0.29869840732194586, -500],  # OHBridge_OTop Dummy (Alpha (110)
    [-599.3929836, 4, 0, 0, 0, 0, -0.8947342143247488, -0.8723293190174063, -603.7210318469869],  # OHBridge_VacTop (Alpha (110))
    [-640.49878946, 4, 0, 4, 0, 0, -0.926725693267018, -0.532342907630269, -644.5146729806254],  # OHBridge_OHTop ****VRange**** (Alpha (110) 
    [-640.49878946, 4, 4, 0, 0, 0, -0.8469360154918735, -0.29869840732194586, -500],  # OHBridge_OTop Dummy (Alpha (110)
    [-640.49878946, 4, 4, 0, 0, 0, -0.8469360154918735, -0.29869840732194586, -500],  # OHBridge_OTop Dummy (Alpha (110)
    [-640.49878946, 4, 4, 0, 0, 0, -0.8469360154918735, -0.29869840732194586, -500],  # OHBridge_OTop Dummy (Alpha (110)
    [-640.49878946, 4, 0, 4, 0, 0, -0.8469360154918735, -0.29869840732194586, -644.6036318649531],  # OHBridge_OHTop (Alpha (110)    
    [-640.49878946, 4, 4, 0, 0, 0, -0.8469360154918735, -0.29869840732194586, -500],  # OHBridge_OTop Dummy (Alpha (110)
]


#### Gibbs Free Energy Corrections for Molecules & Adsorbates ####

### Gas Molecules ###
# 500 eV, O regular, CatHub
E_h2_gas =  -6.77149190
E_h2o_gas = -14.23091949

# (VASP-PBE calculated by Max; T= 300K)
ZPE_H2O = 0.560  # exp. NIST 0.558
ZPE_H2 = 0.268  # exp. NIST 0.273
CV_H2O = 0.103  # From Colin at P = 0.035 bar
CV_H2 = 0.0905
TS_H2O = 0.675  # From Colin at P = 0.035 bar
TS_H2 = 0.408  # at P = 1bar

# Gibbs Free Energy Corrections for Gas Molecules (dE -> dG)
Gcorr_H2O_gas = (ZPE_H2O)+(CV_H2O)-(TS_H2O)
Gcorr_H2_gas = (ZPE_H2)+(CV_H2)-(TS_H2)

### Adsorbates ###
# (VASP-PBE, Top Site on Rutile MnO2 (110), d=0.005; T= 300K)
# Commented Values are from Rutile-MnO2 (110) Bridge Sites d=0.005, OOH data from Alpha-K0.125MnO2(100) Bridge Site d=0.005
ZPE_O = 0.066  # 0.077
ZPE_OH = 0.374  # 0.403
ZPE_OOH = 0.476  # 0.478
CV_O = 0.033  # 0.025
CV_OH = 0.041  # 0.029
CV_OOH = 0.071  # 0.075
TS_O = 0.058  # 0.038
TS_OH = 0.067  # 0.042 
TS_OOH = 0.120  # 0.152

# Gibbs Free Energy Corrections for Adsorbates
Gcorr_O_ads = (ZPE_O)+(CV_O)-(TS_O)
Gcorr_OH_ads = (ZPE_OH)+(CV_OH)-(TS_OH)
Gcorr_OOH_ads = (ZPE_OOH)+(CV_OOH)-(TS_OOH)
Gcorr_H_ads = Gcorr_OH_ads-Gcorr_O_ads
Gcorr_H2O_ads = Gcorr_OH_ads+Gcorr_H_ads


### Function Definitions ####

def addO(pH, E):
    E_O_gas = E_h2o_gas-E_h2_gas
    Gcorr_O_gas = Gcorr_H2O_gas-Gcorr_H2_gas
    return Gcorr_O_ads-(E_O_gas+Gcorr_O_gas)-2*(E+pH*RHE_Shift)

def addOH(pH, E):
    E_OH_gas = E_h2o_gas-0.5*E_h2_gas
    Gcorr_OH_gas = Gcorr_H2O_gas-0.5*Gcorr_H2_gas
    return Gcorr_OH_ads-(E_OH_gas+Gcorr_OH_gas)-(E+pH*RHE_Shift)

def addOOH(pH, E):
    E_OOH_gas = 2*E_h2o_gas-1.5*E_h2_gas
    Gcorr_OOH_gas = 2*Gcorr_H2O_gas-1.5*Gcorr_H2_gas
    return Gcorr_OOH_ads-(E_OOH_gas+Gcorr_OOH_gas)-3*(E+pH*RHE_Shift)

def addH2O(pH, E):
    return Gcorr_H2O_ads-(E_h2o_gas+Gcorr_H2O_gas)

def addH(pH, E):
    E_H_gas = 0.5*E_h2_gas
    Gcorr_H_gas = 0.5*Gcorr_H2_gas
    return Gcorr_H_ads-(E_H_gas+Gcorr_H_gas)+1*(E+pH*RHE_Shift)

## Function to calculate dG (relative to the first entry) for a given coverage (i) at a given pH and voltage (U)
def dg(i, pH, U):
    dg1 = (surfs[i][6]*(U**2) + surfs[i][7]*U + surfs[i][8]) - (surfs[0][6]*(U**2) + surfs[0][7]*U + surfs[0][8]) + surfs[i][1]*addH(pH, U) + surfs[i][2]*addO(pH, U) + surfs[i][3]*addOH(pH, U) + surfs[i][4]*addOOH(pH, U) + surfs[i][5]*addH2O(pH, U)
    return dg1


##### Generate Surface Pourbaix Diagram ######

nsurfs = len(surfs) # Number of Coverages
lowest_surfaces = np.ones((len(E),len(pH)))*-1

# Loop find the coverage with the lowest dG at each value of pH, voltage
pHindex = 0
for i in pH:
    Eindex = 0
    for j in E:
        values = []
        for k in range(nsurfs):
            values.append(dg(k, i, j))
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Eindex][pHindex] = sorted_values[0]
        Eindex+=1
    pHindex+=1

# Defines the colormap used in the surface PBX plot ##
cmapName = 'RdYlBu'

# Creates a meshgrid using each value of pH and voltage
y, x = np.meshgrid(pH,E)

# Plots the lowest energy surface at each value of pH and voltage
plt.pcolormesh(y,x,lowest_surfaces,shading='auto', norm = None, cmap = cmapName, alpha=0.85, vmin=0,vmax=nsurfs)

# Creates plot legend
cmap = plt.get_cmap(cmapName, nsurfs+1)

for k in range(nsurfs): 
    foo = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i H2O-%i)" % (k,
                                          surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4], surfs[k][5])  # Legend Label
    plot([], [], color=cmap(k), linewidth=5, label=foo)

legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1,
         fontsize='x-small', handlelength=3, edgecolor='black')

# Drawing HER and OER Lines on Surface PBX #
pH_forOERLine = np.arange(0, 14, 0.01)
plot(pH_forOERLine, 1.23-pH_forOERLine*RHE_Shift, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 0.95, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$',
        color='blue', rotation=-14, fontsize=10)
plot(pH_forOERLine, 0-pH_forOERLine*RHE_Shift, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, -0.2 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$',
        color='blue', rotation=-14, fontsize=10)

# Saving the Surface PBX Plot #
#fig.savefig(figName+'.pdf', bbox_inches='tight') 
fig.savefig(figName+'.png', bbox_inches='tight') 



#### Plotting Linecut of Surface PBX (dG vs. U) #######

# Define new figure and axes
fig = plt.figure(dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
ax.axis([-1, 2.5, -500., 500.])
ax.tick_params(direction="in")
ax.set_xlabel(r'E (V vs. SHE)')
ax.set_ylabel(r'$\Delta$G (kJ/mol)')

# pH of linecut
pH_linePlot = 14

# Values for X-axis (Voltage vs. SHE)
Es = np.arange(-1.0, 2.7, 0.05)


### Generate Plot ###
for k in range(nsurfs):
    foo = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i H2O-%i)" % (k,
                                          surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4], surfs[k][5])  # Legend Label
    ax.plot(Es, dg(k, pH_linePlot, Es)*kjmol, linestyle ='-', lw=1, c=cmap(k), label=foo)

# Legend
legend(bbox_to_anchor=(0, 1.4), loc=2, borderaxespad=0., ncol=2,
       fontsize='x-small', handlelength=3, edgecolor='black')

# Saving the Linecut Plot #
#fig.savefig(figName+'_Linecut_pH'+str(pH_linePlot)+'.pdf', bbox_inches='tight')
fig.savefig(figName+'_Linecut_pH'+str(pH_linePlot)+'.png', bbox_inches='tight')

# Display the Plot #
# show()
