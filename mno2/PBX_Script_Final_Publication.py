#!/usr/bin/env python

plotname='Pourbaix_PBE+U_U2p75_Mn1E-6_Zn_Epsilon'


import sys, os 
import pymatgen
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion #Decomposes ion into composition object that can handle the charge string
from pymatgen.core import Element #Accesses properties of element inputed (as a string) using an enumeration
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry, MultiEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility #GGA/GGA+U Mixing Scheme
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.plotting import pretty_plot
import sys
import warnings
warnings.filterwarnings('ignore')

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


#This initializes the REST adaptor. Put your own API key in.
mpr = MPRester('b88Lps3sjkIqUfs8GBNDMCL4QeEhUbNW')


##### Constants #####
kJmol = 96.485


print('\n################## Solid Entries ##########################################\n')

pbx_solid_entries = []


##### Delta ######
'''
##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
#'MnO' : -363.4867678050865/kJmol,
#'Mn3O4' : -1294.3840604646095/kJmol, 
#'Mn2O3' : -886.1520794995738/kJmol,

#'MnO2' : -456.4194522641869/kJmol, #Delta
'HMn4O8' : -467.4643620978616*4/kJmol,
'HMn2O4' : -496.0622534546362*2/kJmol, 
#'MnOOH' : -558.8619486482868/kJmol, # (Feitknechtite)
#'Mn(OH)2' : -617.7977867001371/kJmol,  

#'NaMn4O8' : -529.5649780576057*4/kJmol,
#'NaMn2O4' : -595.8992793901739*2/kJmol,
#'NaMnO2' : -700.5005412476619/kJmol, 
'KMn4O8' : -537.5594416232931*4/kJmol,
'KMn2O4' : -587.3354910391488*2/kJmol,
#'LiMn4O8' : -532.3454472589369*4/kJmol, #Delta
#'LiMn2O4' : -613.6138089848864*2/kJmol,
#'LiMnO2' : -751.4300719927867/kJmol,
#'MgMn4O8' : -580.9339849673801*4/kJmol, #Delta
#'MgMn2O4' : -693.5672052188245*2/kJmol,
#'MgMnO2' : -861.4255425333619/kJmol,
#'CaMn4O8' : -610.3130240639866*4/kJmol, #Delta
#'CaMn2O4' : -731.9304954141866*2/kJmol,
#'CaMnO2' : -903.7094643425364/kJmol,
#'SnMn4O8' : -497.73288149911684*4/kJmol, #Delta
#'SnMn2O4' : -522.8814579410965*2/kJmol,
#'SnMnO2' : -472.2519605143566/kJmol,
#'AlMn8O16' : -523.3224635001901*8/kJmol, #Delta
#'AlMn4O8' : -602.475499791193*4/kJmol,
#'AlMn2O4' : -713.4139095869992*2/kJmol,
#'Al2Mn3O6' : -786.4009743780204*3/kJmol,
#'ZnMn4O8' : -506.10477551244895*4/kJmol, #Delta
#'ZnMn2O4' : -545.6534907942117*2/kJmol,
#'ZnMnO2' : -584.1993931540368/kJmol,

## Hydrated ##


#EXP----
'MnO' : -362.898/kJmol, #Dean Exp.
'Mn3O4' : -1283.232/kJmol,
#'Mn2O3' : -881.114/kJmol, #Exp.

'MnO2' : -453.1272/kJmol, #Delta Exp.
'MnOOH' : -543.0832/kJmol, #Delta (Feitknechtite)
'Mn(OH)2' : -615.63376/kJmol, #Delta 
#'MnO2' : -465.138/kJmol, #Beta (rutile) Exp.
#'MnOOH' : -557.7/kJmol, #Beta (manganite)
}
'''


###### Alpha ######
'''
##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
#'MnO' : -363.4867678050865/kJmol,
#'Mn3O4' : -1294.3840604646095/kJmol, 
#'Mn2O3' : -886.1520794995738/kJmol,

'MnO2' : -465.71638408058675/kJmol, #Alpha
'HMn8O16' : -472.454106166424*8/kJmol,
'HMn4O8' : -481.34807889166166*4/kJmol,
'HMn2O4' : -505.35005103608654*2/kJmol, #Alpha
'HMnO2': -531.5152572001871/kJmol, #Alpha
#'H2MnO2': -598.6089087559872/kJmol, #Alpha

'LiMn8O16' : -501.8760566661115*8/kJmol,
'LiMn4O8' : -534.046884975637*4/kJmol,
'LiMn2O4' : -607.1731284125865*2/kJmol, #Alpha
#'NaMn8O16' : -504.10445822069596*8/kJmol,
#'NaMn4O8' : -526.4786062291054*4/kJmol,
#'NaMn2O4' : -469.76339299277396*2/kJmol, #Alpha
#'KMn8O16' : -509.4734770692398*8/kJmol,
#'KMn4O8' : -533.9320147196931*4/kJmol,
#'MgMn8O16' : -517.3045518816087*8/kJmol,
#'MgMn4O8' : -559.2256135152303*4/kJmol,
#'MgMn2O4' : -660.1020930984243*2/kJmol, #Alpha
#'CaMn8O16' : -536.5552922241368*8/kJmol,
#'CaMn4O8' : -582.7519525943867*4/kJmol,
#'CaMn2O4' : -702.2623507448367*2/kJmol, #Alpha
#'ZnMn8O16' : -481.6435639588431*8/kJmol,
#'ZnMn4O8' : -491.943349802549*4/kJmol,
#'ZnMn2O4' : -526.8540485016116*2/kJmol, #Alpha
#'SnMn8O16' : -491.5767605229518*8/kJmol,
#'SnMn4O8' : -501.0979261288666*4/kJmol,
#'SnMn2O4' : -549.0506011012466*2/kJmol, #Alpha
#'AlMn8O16' : -515.6442499154402*8/kJmol,
#'AlMn4O8' : -561.6157473604427*4/kJmol, #Alpha

## Hydrated ##
#'H3Mn8O17' : -472.59636900228037*8/kJmol,
#'H3Mn4O9' : -476.95021838857457*4/kJmol,
'H2LiMn8O17' : -500.1834951321678*8/kJmol,
'H2LiMn4O9' : -543.3213011427991*4/kJmol,
#'H2NaMn8O17' : -498.90326886910225*8/kJmol,
#'H2KMn8O17' : -500.91306016854605*8/kJmol,
#'H2MgMn8O17' : -517.964338065915*8/kJmol,
#'H2MgMn4O9' : -577.3186978968427*4/kJmol,
#'H2CaMn8O17' : -537.453762036593*8/kJmol,
#'H2ZnMn8O17' : -481.9316580830989*8/kJmol,
#'H2ZnMn4O9' : -514.3705809078614*4/kJmol,


#EXP----
#'MnO' : -362.898/kJmol, #Dean Exp.
'Mn3O4' : -1283.232/kJmol,
'Mn2O3' : -881.114/kJmol, #Exp.

#'MnO2' : -465.138/kJmol, #Beta (rutile) Exp.
#'MnOOH' : -557.7/kJmol, #Beta (manganite)
'Mn(OH)2' : -615.63376/kJmol, #Delta 
}
'''


##### Lowest E #####

##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
"Zn": 0,
# 'MnO2' : -465.138/kJmol, #Beta(rutile)Exp
'MnO2' : -440.16/kJmol, #Epsilon PBE+U=2.75
'Mn2O3' : -881.114/kJmol, #Exp.
'Mn3O4' : -1283.232/kJmol, #Exp.
'MnOOH' : -557.7/kJmol, #Beta (manganite) Exp
'Mn(OH)2' : -615.63376/kJmol, #Delta Exp

# 'HMn8O16' : -472.454106166424*8/kJmol,#Alpha
# 'HMn4O8' : -482.1357496267614*4/kJmol, #Rams
# 'HMn2O4' : -505.35005103608654*2/kJmol, #Alpha
#'MnOOH' : -555.836682699987/kJmol, #Beta (manganite)
#'MnOOH' : -558.8619486482868/kJmol, # (Feitknechtite)
#'Mn(OH)2' : -617.7977867001371/kJmol #Delta
#'LiMn8O16' : -501.8760566661115*8/kJmol, #Alpha
#'LiMn4O8' : -541.4949017528368*4/kJmol, #Spinel
#'LiMn2O4' : -628.4736761018866*2/kJmol, #Spinel
#'LiMnO2' : -751.4300719927867/kJmol, #Delta
#'NaMn8O16' : -504.10445822069596*8/kJmol,#Alpha
#'NaMn4O8' : -529.5649780576057*4/kJmol,#Delta
#'NaMn2O4' : -595.8992793901739*2/kJmol,#Delta
#'NaMnO2' : -700.5005412476619/kJmol, #Delta
# 'KMn8O16' : -509.4734770692398*8/kJmol, #Alpha
# 'KMn4O8' : -537.5594416232931*4/kJmol, #Delta
# 'KMn2O4' : -587.3354910391488*2/kJmol,
#'MgMn8O16' : -517.3045518816087*8/kJmol, #Alpha
#'MgMn4O8' : -580.9339849673801*4/kJmol, #Delta
#'MgMn2O4' : -709.6991141733745*2/kJmol, #Spinel
#'MgMnO2' : -861.4255425333619/kJmol, #Delta
#'CaMn8O16' : -536.5552922241368*8/kJmol,#Alpha
#'CaMn4O8' : -610.3130240639866*4/kJmol, #Delta
#'CaMn2O4' : -731.9304954141866*2/kJmol, #Delta
#'CaMnO2' : -903.7094643425364/kJmol, #Delta
'ZnMn8O16' : -481.6435639588431*8/kJmol, #Alpha
'ZnMn4O8' : -516.324975188399*4/kJmol, #Spinel
'ZnMn2O4' : -584.5695268490117*2/kJmol, #Spinel
'ZnMnO2' : -584.1993931540368/kJmol, #Delta
#'AlMn8O16' : -525.3090491264901*8/kJmol, #Spinel
#'AlMn4O8' : -602.475499791193*4/kJmol, #Delta
#'AlMn2O4' : -737.8049654167493*2/kJmol, #Spinel
#'Al2Mn3O6' : -786.4009743780204*3/kJmol, #Delta

#Precipitated Ion Solids
#'LiOH' : -438.954/kJmol,
#'NaOH' : -379.737/kJmol,
# 'KOH' : -378.858/kJmol,
#'Mg(OH)2' : -833.644/kJmol,
#'Ca(OH)2' : -898.470/kJmol,
'ZnO' : -320.476/kJmol,
#'Al(OH)3' : -1138.706/kJmol,
}

'''
##### Beta #####

##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
#'MnO' : -363.4867678050865/kJmol,
#'Mn3O4' : -1294.3840604646095/kJmol, 
#'Mn2O3' : -886.1520794995738/kJmol,

#'MnO2' : -458.58976848688667/kJmol, #Beta (rutile) 
'HMn4O8' : -472.5964185448616*4/kJmol, #Beta
'HMn2O4' : -493.25367930843635*2/kJmol, #Beta
'MnOOH' : -555.836682699987/kJmol, #Beta (manganite)
#'Mn(OH)2' : -599.8865061998371/kJmol, #Beta  

#'LiMn8O16' : -484.5598995052118*8/kJmol,
#'LiMn4O8' : -506.56561817923676*4/kJmol,
#'LiMn2O4' : -574.4246200356864*2/kJmol, #Beta
#'LiMnO2' : -728.8622787352867/kJmol, #Beta
#'MgMn8O16' : -495.7297272476587*8/kJmol,
#'MgMn4O8' : -531.6189905122803*4/kJmol,
'MgMn2O4' : -665.0217503609745*2/kJmol, #Beta
#'ZnMn8O16' : -461.84926386639313*8/kJmol,
#'ZnMn4O8' : -465.018793737349*4/kJmol,
#'ZnMn2O4' : -543.6049193445118*2/kJmol, #Beta
#'SnMn8O16' : -472.334273583102*8/kJmol,
#'SnMn4O8' : -499.3569314918664*4/kJmol,
#'SnMn2O4' : -468.56482774319664*2/kJmol, #Beta
'AlMn8O16' : -519.11347132939*8/kJmol,
'AlMn4O8' : -575.8952658860926*4/kJmol,
'AlMn2O4' : -712.9212156884494*2/kJmol, #Beta

#EXP----
'MnO' : -362.898/kJmol, #Dean Exp.
'Mn3O4' : -1283.232/kJmol,
#'Mn2O3' : -881.114/kJmol, #Exp.

'MnO2' : -465.138/kJmol, #Beta (rutile) Exp.
#'MnOOH' : -557.7/kJmol, #Beta (manganite)
'Mn(OH)2' : -615.63376/kJmol, #Delta 
}
'''
'''
##### Spinel #####

##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
#'MnO' : -363.4867678050865/kJmol,
#'Mn3O4' : -1294.3840604646095/kJmol, 
#'Mn2O3' : -886.1520794995738/kJmol,

'MnO2' : -452.5362745007868/kJmol, #Spinel
'HMn4O8' : -472.31159289516137*4/kJmol, #Spinel
'HMn2O4' : -505.0423256364865*2/kJmol, #Spinel
'MnOOH' : -531.0296617030871/kJmol, #Spinel
'Mn(OH)2' : -608.5305288454871/kJmol, #Spinel 

'LiMn4O8' : -541.4949017528368*4/kJmol, #Spinel
'LiMn2O4' : -628.4736761018866*2/kJmol, #Spinel
'LiMnO2' : -750.4065822691867/kJmol, #Spinel
#'NaMn4O8' : -519.1123614542555*4/kJmol, #Spinel
#'NaMn2O4' : -588.1905234786739*2/kJmol, #Spinel
#'NaMnO2' : -691.859036860512/kJmol, #Spinel
#'KMn4O8' : -474.500073046543*4/kJmol, #Spinel
#'KMn2O4' : -524.267288295799*2/kJmol, #Spinel
#'KMnO2' : -585.7509949541117/kJmol, #Spinel
#'MgMn4O8' : -580.1749278238802*4/kJmol, #Spinel
#'MgMn2O4' : -709.6991141733745*2/kJmol, #Spinel
#'MgMnO2' : -860.2383283781118/kJmol, #Spinel
#'CaMn4O8' : -588.4935578531365*4/kJmol, #Spinel
#'CaMn2O4' : -731.5925759986868*2/kJmol, #Spinel
#'CaMnO2' : -903.3022117708865/kJmol, #Spinel
#'ZnMn4O8' : -516.324975188399*4/kJmol, #Spinel
#'ZnMn2O4' : -584.5695268490117*2/kJmol, #Spinel
#'ZnMnO2' : -582.5095374967365/kJmol, #Spinel
#'AlMn8O16' : -525.3090491264901*8/kJmol, #Spinel
#'AlMn4O8' : -597.5234033067925*4/kJmol, #Spinel
#'AlMn2O4' : -737.8049654167493*2/kJmol, #Spinel
#'Al2Mn3O6' : -781.0474925748204*3/kJmol

#EXP----
#'MnO' : -362.898/kJmol, #Dean Exp.
'Mn3O4' : -1283.232/kJmol,
#'Mn2O3' : -881.114/kJmol, #Exp.

#'MnO2' : -465.138/kJmol, #Beta (rutile) Exp.
#'MnOOH' : -557.7/kJmol, #Beta (manganite)
#'Mn(OH)2' : -615.63376/kJmol, #Delta 
}
'''

##### Ramsdellite #####
'''
##PBE Energies U = 2.75
ion_dict_solids_expt = {
'Mn': 0,
#'MnO' : -363.4867678050865/kJmol,
#'Mn3O4' : -1294.3840604646095/kJmol, 
#'Mn2O3' : -886.1520794995738/kJmol,

#'MnO2' : -463.79371590738657/kJmol, #Rams
'HMn4O8' : -482.1357496267614*4/kJmol, 
'HMn2O4' : -490.2444974096364*2/kJmol, #Rams
#'MnOOH' : -557.9869107456868/kJmol, #Rams
#'Mn(OH)2' : -607.1812353278372/kJmol, #Rams 

#'LiMn4O8' : -536.0279005228867*4/kJmol, #Rams
#'LiMn2O4' : -611.6320765540867*2/kJmol, #Rams
#'LiMnO2' : -737.5131344486368/kJmol, #Rams
#'NaMn4O8' : -512.6574290826056*4/kJmol, #Rams
#'NaMn2O4' : -567.2773727128739*2/kJmol, #Rams
#'NaMnO2' : -663.8578114342619/kJmol, #Rams
#'MgMn4O8' : -572.41342259803*4/kJmol, #Rams
#'MgMn2O4' : -676.2962011082244*2/kJmol, #Rams
#'MgMnO2' : -835.7865593621118/kJmol, #Rams
#'CaMn4O8' : -553.0212190703867*4/kJmol, #Rams
#'CaMn2O4' : -694.0966370490868*2/kJmol, #Rams
#'CaMnO2' : -899.5171303069865/kJmol, #Rams
#'ZnMn4O8' : -503.0294425327992*4/kJmol, #Rams
#'ZnMn2O4' : -549.5152286257618*2/kJmol, #Rams
#'ZnMnO2' : -569.4497553252867/kJmol, #Rams
'AlMn8O16' : -523.60034512444*8/kJmol, #Rams
'AlMn4O8' : -593.0169613888927*4/kJmol, #Rams
'AlMn2O4' : -716.6765701238493*2/kJmol, #Rams
'Al2Mn3O6' : -800.1189251175703*3/kJmol #Rams

#EXP----
'MnO' : -362.898/kJmol, #Dean Exp.
'Mn3O4' : -1283.232/kJmol,
'Mn2O3' : -881.114/kJmol, #Exp.

#'MnO2' : -465.138/kJmol, #Beta (rutile) Exp.
'MnO2' : -462.06666499999994/kJmol, #Rams
'MnOOH' : -557.7/kJmol, #Beta (manganite)
'Mn(OH)2' : -615.63376/kJmol, #Delta 
}
'''

conc_dict = {}

for key in ion_dict_solids_expt:
	comp = Ion.from_formula(key)
	energy = ion_dict_solids_expt[key]
	pbx_entry_ion = PourbaixEntry(ComputedEntry(comp, energy))
	pbx_solid_entries.append(pbx_entry_ion)
	conc_dict[key] = 1

for a in pbx_solid_entries:
	print(a, a.concentration)
	#print(a)




print('\n################## Ion Entries ##########################################\n')

pbx_ion_entries = []

# Helper function to convert new API format to old format
def convert_ion_data(new_format_list):
    """Convert new API format to old format expected by the code"""
    old_format = []
    for item in new_format_list:
        if isinstance(item, dict):
            # New format: {'identifier': ..., 'formula': ..., 'data': {...}}
            if 'identifier' in item or 'formula' in item:
                name = item.get('identifier') or item.get('formula')
                data = item.get('data', {})
                # Energy is in kJ/mol, convert to eV
                energy_kjmol = data.get('ΔGᶠ', {}).get('value', 0)
                energy_ev = energy_kjmol / kJmol if energy_kjmol else 0
                
                old_format.append({
                    'Name': name,
                    'Energy': energy_ev,
                    'Reference Solid': data.get('RefSolid', ''),
                    'Reference solid energy': data.get('ΔGᶠRefSolid', {}).get('value', 0) / kJmol if isinstance(data.get('ΔGᶠRefSolid'), dict) else 0,
                    'Source': data.get('reference', 'MP'),
                    'Major_Elements': [data.get('MajElements', '')]
                })
            # Old format: {'Name': ..., 'Energy': ...}
            elif 'Name' in item and 'Energy' in item:
                old_format.append(item)
    return old_format

# Use get_ion_reference_data_for_chemsys (works with new API)
try:
    ion_dict_Mn_raw = mpr.get_ion_reference_data_for_chemsys(["Mn"])
    if not isinstance(ion_dict_Mn_raw, list):
        ion_dict_Mn_raw = list(ion_dict_Mn_raw) if hasattr(ion_dict_Mn_raw, '__iter__') else []
    
    # Convert to old format
    ion_dict_Mn = convert_ion_data(ion_dict_Mn_raw)
    
    for i,dummy in enumerate(ion_dict_Mn):
        print(ion_dict_Mn[i])
    print()
except Exception as e:
    print(f"Error getting Mn ion reference data: {e}")
    ion_dict_Mn = []

# Get K reference data
try:
    ion_dict_X_raw = mpr.get_ion_reference_data_for_chemsys(["Zn"])
    if not isinstance(ion_dict_X_raw, list):
        ion_dict_X_raw = list(ion_dict_X_raw) if hasattr(ion_dict_X_raw, '__iter__') else []
    
    # Convert to old format
    ion_dict_X = convert_ion_data(ion_dict_X_raw)
    
    for i,dummy in enumerate(ion_dict_X):
        print(ion_dict_X[i])
    print()
except Exception as e:
    print(f"Error getting K ion reference data: {e}")
    ion_dict_X = []


ion_dict = ion_dict_Mn + ion_dict_X


for id in ion_dict_X:
    comp = Ion.from_formula(id['Name']) # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
    energy = id['Energy'] #+ ion_correction * factor
    conc_dict[id['Name']] = 1
    print(id['Name'], comp, energy)
    pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy, id['Name']))
    #if pbx_entry_ion.name not in ["HZnO2[-]","ZnO2[2-]","Zn(OH)2(aq)","ZnOH[+]"]:
    #if pbx_entry_ion.name not in ["HZnO2[-]","ZnO2[2-]","Zn(OH)2(aq)"]:
    #if pbx_entry_ion.name not in ["HZnO2[-]","ZnO2[2-]"]:
    #if pbx_entry_ion.name not in ["AlO2[-]","AlOH[2+]","Al(OH)4[-]"]:
    #if pbx_entry_ion.name not in ["LiOH(aq)"]:
    #if pbx_entry_ion.name not in ["MgOH[+]"]:
    pbx_ion_entries.append(pbx_entry_ion)

#pbx_ion_entries.append(PourbaixEntry(IonEntry(Ion.from_formula('Zn(OH)4[2-]'), -8.94036, 'Zn(OH)4[2-]')))
#conc_dict['Zn(OH)4[2-]'] = 1


for id in ion_dict_Mn:
    comp = Ion.from_formula(id['Name']) # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
    energy = id['Energy'] #+ ion_correction * factor
    conc_dict[id['Name']] = 1e-6
    print(id['Name'], comp, energy)
    pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy, id['Name']))
    if pbx_entry_ion.name not in ["HCoO3[2-]"]:
        pbx_ion_entries.append(pbx_entry_ion)

conc_dict["Mn2+"]=1
conc_dict["Zn2+"]=1

print()
all_entries = pbx_solid_entries + pbx_ion_entries

conc_dict["Zn"]=1

for a in conc_dict:
    print(a, conc_dict[a])

def replot(all_entries):
    pourbaix = PourbaixDiagram(all_entries, comp_dict = {"Mn": 0.5, "Zn": 0.5}, conc_dict=conc_dict)
    plotter = PourbaixPlotter(pourbaix)
    ax = plotter.get_pourbaix_plot(limits=[[0, 16],[-1.5, 2.0]],label_domains=True, show_neutral_axes=False)
    f = ax.figure
    f.set_size_inches((12*1.5, 7*1.5))  # 또는 cm 단위: (30/2.54, 20/2.54)
    
    # xlabel과 ylabel의 fontsize 변경
    ax.set_xlabel(ax.get_xlabel(), fontsize=28)  # 원하는 fontsize로 변경 가능
    ax.set_ylabel(ax.get_ylabel(), fontsize=28)  # 원하는 fontsize로 변경 가능
    
    # xtick과 ytick의 fontsize 변경
    ax.tick_params(axis='x', labelsize=24)  # x축 tick label fontsize
    ax.tick_params(axis='y', labelsize=24)  # y축 tick label fontsize
    
    # entry를 구분하는 선의 두께 변경
    for line in ax.lines:
        line.set_linewidth(1.0)  # 원하는 선 두께로 변경 가능 (예: 0.5, 1.0, 1.5, 2.0 등)
    
    # entry 이름(도메인 레이블)의 fontsize 변경
    for text in ax.texts:
        text.set_fontsize(24)  # 원하는 fontsize로 변경 가능
    
    plt.tight_layout()

    # plt.savefig(plotname+'.pdf')
    plt.savefig(plotname+'.pdf', dpi=300, bbox_inches='tight', transparent=True)
    plt.show()

    # plt2=plotter.get_pourbaix_plot_colorfill_by_element(limits=[[-2, 16],[-3, 3]])
    # plt2.savefig(plotname+'colorfill.pdf')
    # plt2.savefig(plotname+'colorfill.png')
    # plt2.show()
    
    #plotter.show(filename=plotname+'convex_hull_plot.eps')
    #plotter.show()


replot(all_entries)

'''
#### Constants ####
pH = 14
R = 8.314
T = 298.15
muH2O = -237.141/kJmol #Barin dG H2O(l)
mupH = -1*R*T/(kJmol*1000)*np.log(10)*pH

#def pHcut(pH, V1, V2):
'''






