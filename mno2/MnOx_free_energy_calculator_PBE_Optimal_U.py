#!/usr/bin/python

kjmol=96.485 #eV/atom to kj/mol
T=298.15
convert_dS=T/(1000*kjmol) #J/(mol*K) to eV/atom


############################################################################################


print('\n############# Enthalpy and Free Energy of Molecular References ########################################')
print('All computational values from PBE, afm Mn ordering')


#PBE 600eV
h2o=-14.23919983
h2=-6.77409008
o2=-9.89038823


#Check how SCAN matches experiment (PBE very close to exp.)
h2ozpe=0.5584   #zpeh2o=0.560   #exp. NIST 0.558
h2ocp=0.10260 	#NIST 0.10      #cvh2o=0.103

h2zpe=0.27283   #zpeh2=0.268  	#exp. NIST 0.273
h2cp=0.087785   #=NIST  0.09    #cvh2=0.0905

o2zpe=0.098
o2cp=0.089962


#Enthalpy(dH) of H2O and H2 from PBE at 300K
h_ref=(h2+h2zpe+h2cp)/2
h2o_ref=h2o+h2ozpe+h2ocp

#Energy(dE) & Enthalpy(dH) of O2 from PBE at 300K
o_ref1_de=(o2)/2 #dE (Energy)
o_ref1=(o2+o2zpe+o2cp)/2 #dH (Enthalpy)

#Energy(dE) & Enthalpy(dH) of O2 from H2O/H2 PBE values at 300K
#2.506 = -241.81*kjmol (96.485)
o_ref2_de=h2o_ref-2*h_ref +2.506 -(o2zpe+o2cp)/2  #dE (Energy) relative to -241.81/kjmol
o_ref2=h2o_ref-2*h_ref +2.506  #dH fitting relative to -241.81/kjmol

print('\nError in dH of O2 PBE value: ',h2o_ref-2*h_ref-o_ref1 -(-241.81/kjmol), ' (if not close to zero: error in O2 PBE energy)')
print('Error in dH of O2 from H2O/H2 PBE values: ',h2o_ref-2*h_ref-o_ref2 -(-241.81/kjmol), ' (should be close to zero!)')

#Standard Entropy of O2/H2 gas (Experimental Values)
S_o2_gas_exp=205.2 #exp.
S_h2_gas_exp=130.68 #exp.

#T*dS of O2/H2O/H2
tsh2o=0.675 	#From Colin at P = 0.035 bar 
tsh2_exp=S_h2_gas_exp*convert_dS 	#at P = 1bar
tso2_exp=S_o2_gas_exp*convert_dS


#check H2O thermochemistry
print('\ndH of molecular references at 300K:')
print('1/2 O2 derived from H2O/H2: ',o_ref2)
print('1/2 O2 derived from PBE: ', o_ref1)
print('1/2 H2 derived from PBE: ', h_ref)
print('H2O derived from PBE: ', h2o_ref)


print('\nError in dG of O2 PBE value:     ', (h2o_ref -2*h_ref -o_ref1 -(tsh2o-tsh2_exp-0.5*tso2_exp)) -(-237.140/kjmol), ' (if not close to zero: error in O2 PBE energy)')
print('Error in dG of O2 from H2O/H2 PBE values: ', (h2o_ref -2*h_ref -o_ref2 -(tsh2o-tsh2_exp-0.5*tso2_exp)) -(-237.140/kjmol), ' (should be close to zero!)')





############################################################################################

print('\n############# Enthalpy and Free Energy of MnOxHy System ########################################')

#Units

kcalTokJ = 4.184


################################
#Experimental Enthalpies
dh_B_MnO2_exp=-520.029/kjmol 		#Barin 
dh_R_MnO2_exp=5.4+dh_B_MnO2_exp*kjmol 	#Fritsch & Navrotsky, 1998  #+/- 3.1kJ/mol 
dh_Mn2O3_exp=-959.002 			#Barin
dh_Mn3O4_exp=-1387.799/kjmol 		#Barin
dh_MnO_exp=-385.221 			#Barin
dh_D_MnOH2_exp=-163.2*kcalTokJ/kjmol 	#Oswald & Asper (from Nb1 originally)
dh_B_MnOOH_exp=(-12.6+0.5*dh_Mn2O3_exp-0.5*285.830)/kjmol #Fritsch & Navrotsky, 1996 	#H2O from Barin #+/- 0.9kJ/mol
dh_R_MnOOH_exp=-10.9+0.5*dh_Mn2O3_exp-0.5*285.830 	  #Fritsch & Navrotsky, 1996 	#H2O from Barin #+/- 1.7kJ/mol



################################
#Experimental Free Energies
dg_B_MnO2_exp=-465.138 		#Pyrolusite 	#Barin (kJ/mol) #-111.17*kcalTokJ=-465.135 Hem & Lind, 1983
dg_D_MnO2_exp=-108.3*kcalTokJ 	#Birnessite 	#Hem & Lind, 1983 #probably intercalated ions/water
dg_R_MnO2_exp=-4.789*kjmol	#Ramsdellite	#Chen & Schelhas, 2018
dg_G_MnO2_exp=-109.1*kcalTokJ 	#EMD 	  	#Hem & Lind, 1983
dg_B_MnOOH_exp=-133.3*kcalTokJ 	#Manganite 	#Hem & Lind, 1983
dg_D_MnOOH_exp=-129.8*kcalTokJ 	#Feitknechtite 	#Hem & Lind, 1983
dg_Mn2O3_exp=-881.114 		#Bixbyite	#Barin
dg_Mn3O4_exp=-1283.232 		#Hausmannite	#Barin (kJ/mol) #-306.7*kcalTokJ=-1283.233 Hem & Lind, 1983
dg_D_MnOH2_exp=-147.14*kcalTokJ #Pyrochroite	#Hem & Lind, 1983
dg_MnO_exp=-362.898 		#Manganosite	#Barin



################################
#Experimental Entropies
S_MnO2_solid_exp = 53.049 	#Barin 		#52.75 Robie and Hemingway (1995)
S_Mn2O3_exp = 110.499 		#Barin J/molK 	#113.7 Robie and Hemingway (1995)
S_Mn3O4_exp = 155.599 		#Barin
S_D_MnOH2_exp = 23.7*kcalTokJ 	#Oswald & Asper
S_MnO_exp = 59.710 		#Barin 
S_Mn_metal_exp = 32.008  	#Barin
dS_form_B_MnOOH_exp = (dg_B_MnOOH_exp-dh_B_MnOOH_exp*kjmol)/-298.15*1000 	#-185.43283582089546J/molK
S_B_MnOOH_exp = dS_form_B_MnOOH_exp+S_Mn_metal_exp+S_o2_gas_exp+S_h2_gas_exp/2 #Calculating S_B_MnOOH_exp from dH and dG experimental values
print(S_B_MnOOH_exp)	


################################
#Zero Point Energies (eV/fu) (PBE+U, U=3.75)

MnO2_ZPE = 0.2867
R_MnO2_ZPE = 0.2996
D_MnO2_ZPE = 0.2988
MnOOH_ZPE = 0.5779
MnOH2_ZPE = 0.8617
Mn3O4_ZPE = 0.6573
Mn2O3_ZPE = 0.4937


S_MnO2_PBE = 0.0006337 #eV/atom/K close to 53.049/1000/96.485=0.00054982
S_MnOOH_PBE = 0.0005301 #eV/atom/K 
S_MnOH2_PBE = 0.0008545 #eV/atom/K close to 23.7*4.184/1000/96.485=0.00102773


################################
#MnOxHy Bulk Energies (eV/fu) (PBE+U, U=3.75)
Mn_metal = -5.929491398 #PBE
B_MnOH2 = -30.2613637+MnOH2_ZPE


'''
#####
#U=3.75
B_MnO2 = -21.20169628+MnO2_ZPE
R_MnO2 = -21.27653229+R_MnO2_ZPE 
B_MnOOH = -26.18287205+MnOOH_ZPE
D_MnOH2 = -30.45498958+MnOH2_ZPE #-30.4051931
Mn3O4 = -53.52859596+Mn3O4_ZPE
Mn2O3 = -37.50301343+Mn2O3_ZPE
D_MnO2 = -21.21466879+D_MnO2_ZPE
D_MnOOH = -26.20729815+MnOOH_ZPE
'''

'''
#####
#U=3.5
B_MnO2 = -21.37084707+MnO2_ZPE
R_MnO2 = -21.44395586+R_MnO2_ZPE 
B_MnOOH = -26.29158087+MnOOH_ZPE
D_MnOH2 = -30.50065703+MnOH2_ZPE #-30.45180171
Mn3O4 = -53.79881115+Mn3O4_ZPE
Mn2O3 = -37.71211329+Mn2O3_ZPE
D_MnO2 = -21.38066379+D_MnO2_ZPE
D_MnOOH = -26.32216084+MnOOH_ZPE
'''

'''
#####
#U=3.25
B_MnO2 = -21.53753521+MnO2_ZPE
R_MnO2 = -21.60843412+R_MnO2_ZPE 
B_MnOOH = -26.39744157+MnOOH_ZPE
D_MnOH2 = -30.54745695+MnOH2_ZPE #-30.49898729
Mn3O4 = -54.07480468+Mn3O4_ZPE
Mn2O3 = -37.92552693+Mn2O3_ZPE
D_MnO2 = -21.54076578+D_MnO2_ZPE
D_MnOOH = -26.42827484+MnOOH_ZPE
'''

'''
#####
#U=3
B_MnO2 = -21.70537094+MnO2_ZPE
R_MnO2 = -21.77482841+R_MnO2_ZPE 
B_MnOOH = -26.50491918+MnOOH_ZPE
D_MnOH2 = -30.59540937+MnOH2_ZPE #-30.54726124
Mn3O4 = -54.35544474+Mn3O4_ZPE
Mn2O3 = -38.14158369+Mn2O3_ZPE
D_MnO2 = -21.70249795+D_MnO2_ZPE
D_MnOOH = -26.53596001+MnOOH_ZPE
'''

'''
#####
#U=2.93
B_MnO2 = -21.75350094+MnO2_ZPE
R_MnO2 = -21.82172273+R_MnO2_ZPE 
B_MnOOH = -26.53530736+MnOOH_ZPE
D_MnOH2 = -30.60902279+MnOH2_ZPE #-30.56099131
Mn3O4 = -54.43524617+Mn3O4_ZPE
Mn2O3 = -38.20380821+Mn2O3_ZPE
D_MnO2 = -21.74807847+D_MnO2_ZPE
D_MnOOH = -26.56634488+MnOOH_ZPE
'''

#####
#U=2.75
B_MnO2 = -21.87625236+MnO2_ZPE
R_MnO2 = -21.94308766+R_MnO2_ZPE 
B_MnOOH = -26.61405836+MnOOH_ZPE
D_MnOH2 = -30.64216489+MnOH2_ZPE #-30.59671554
Mn3O4 = -54.64187765+Mn3O4_ZPE
Mn2O3 = -38.36185625+Mn2O3_ZPE
D_MnO2 = -21.86585854+D_MnO2_ZPE
D_MnOOH = -26.64541314+MnOOH_ZPE
E_MnO2 = -21.6852213333+MnO2_ZPE

'''
#####
#U=2.5
B_MnO2 = -22.04827663+MnO2_ZPE
R_MnO2 = -22.11297896+R_MnO2_ZPE 
B_MnOOH = -26.72493741+MnOOH_ZPE
D_MnOH2 = -30.69474626+MnOH2_ZPE  #-30.64738154
Mn3O4 = -54.93201978+Mn3O4_ZPE
Mn2O3 = -38.58684185+Mn2O3_ZPE
D_MnO2 = -22.03034906+D_MnO2_ZPE
D_MnOOH = -26.75653534+MnOOH_ZPE
'''

'''
#####
#U=2.25
B_MnO2 = -22.22136572+MnO2_ZPE
R_MnO2 = -22.28410511+R_MnO2_ZPE 
B_MnOOH = -26.8367996+MnOOH_ZPE
D_MnOH2 = -30.74627147+MnOH2_ZPE
Mn3O4 = -55.23016574+Mn3O4_ZPE
Mn2O3 = -38.81345113+Mn2O3_ZPE
D_MnO2 = -22.19682548+D_MnO2_ZPE
D_MnOOH = -26.86949954+MnOOH_ZPE
'''

'''
#####
#U=0
B_MnO2 = -23.87915056+MnO2_ZPE
R_MnO2 = -23.92124592+R_MnO2_ZPE 
B_MnOOH = -27.94368017+MnOOH_ZPE
D_MnOH2 = -31.2693512+MnOH2_ZPE #-31.22677161
Mn3O4 = -58.22526961+Mn3O4_ZPE
Mn2O3 = -41.06658923+Mn2O3_ZPE
D_MnO2 = -23.7645086+D_MnO2_ZPE
D_MnOOH = -27.9768031+MnOOH_ZPE
'''

################################
#Schemes to produce better fit between dG of experiment and dG of computation
Mn_metal_fit_Mn4=B_MnO2-(2*o_ref2)-dh_B_MnO2_exp  #fit to exp dH of MnO2
Mn_metal_fit_Mn4R=R_MnO2-(2*o_ref2)-dh_R_MnO2_exp/kjmol  #fit to exp dH of R-MnO2
Mn_metal_fit_Mn3_Mn2O3=(Mn2O3-(3*o_ref2)-dh_Mn2O3_exp/kjmol)/2  #fit to exp dH of Mn2O3
Mn_metal_fit_Mn3_MnOOH=B_MnOOH-(2*o_ref2+h_ref)-dh_B_MnOOH_exp #fit to exp dH of MnOOH
Mn_metal_fit_Mn2=D_MnOH2-(2*o_ref2+2*h_ref)-dh_D_MnOH2_exp #fit to exp dH of MnOH2
Mn_metal_fit_Mn23=(Mn3O4-(4*o_ref2)-dh_Mn3O4_exp)/3 #fit to exp dH of Mn3O4

print('Calculated Mn metal ref: ', Mn_metal)
print('Mn ref fit to B_MnO2 dH : ', Mn_metal_fit_Mn4)
print('Mn ref fit to R_MnO2 dH : ', Mn_metal_fit_Mn4R)
print('Mn ref fit to B_MnOOH dH: ', Mn_metal_fit_Mn3_MnOOH)
print('Mn ref fit to B_Mn2O3 dH: ', Mn_metal_fit_Mn3_Mn2O3)
print('Mn ref fit to D_MnOH2 dH : ',  Mn_metal_fit_Mn2)
print('Mn ref fit to Mn3O4 dH : ', Mn_metal_fit_Mn23)

Mn_metal_fit_avg = (Mn_metal_fit_Mn4+Mn_metal_fit_Mn3_MnOOH+Mn_metal_fit_Mn3_Mn2O3+Mn_metal_fit_Mn2+Mn_metal_fit_Mn23)/5
print(Mn_metal_fit_avg)

Mn_metal_fit_avgR = (Mn_metal_fit_Mn4+Mn_metal_fit_Mn4R+Mn_metal_fit_Mn3_MnOOH+Mn_metal_fit_Mn3_Mn2O3+Mn_metal_fit_Mn2+Mn_metal_fit_Mn23)/6


#Using perfect energy for metal relative to experimental dH of MnO2 
Mn_metal4 = Mn_metal_fit_avg
Mn_metal2 = Mn_metal_fit_avg
Mn_metal3 = Mn_metal_fit_avg
Mn_metal23 = Mn_metal_fit_avg

#Mn_metal4 = Mn_metal
#Mn_metal2 = Mn_metal
#Mn_metal3 = Mn_metal
#Mn_metal23 = Mn_metal




################################
print('\nCalculate dHs now')
dh_B_MnO2=B_MnO2-(2*o_ref2+Mn_metal4)
print('B-MnO2: dH=' ,dh_B_MnO2, 'eV  ',dh_B_MnO2*kjmol, 'kj/mol   exp: ', dh_B_MnO2_exp*kjmol, 'kjmol')

dh_R_MnO2=R_MnO2-(2*o_ref2+Mn_metal4)
print('R-MnO2: dH=' ,dh_R_MnO2, 'eV  ',dh_R_MnO2*kjmol, 'kj/mol   exp: ', dh_R_MnO2_exp, 'kjmol')

dh_D_MnO2=D_MnO2-(2*o_ref2+Mn_metal4)
print('D-MnO2: dH=' ,dh_D_MnO2, 'eV  ', dh_D_MnO2*kjmol, 'kj/mol')

dh_E_MnO2=E_MnO2-(2*o_ref2+Mn_metal4)
print('E-MnO2: dH=' ,dh_E_MnO2, 'eV  ', dh_E_MnO2*kjmol, 'kj/mol')

dh_B_MnOOH=B_MnOOH-(2*o_ref2+Mn_metal3+h_ref)
print('B-MnOOH: dH=',dh_B_MnOOH, 'eV  ',dh_B_MnOOH*kjmol, 'kj/mol   exp: ', dh_B_MnOOH_exp*kjmol, 'kjmol')

dh_D_MnOOH=D_MnOOH-(2*o_ref2+Mn_metal3+h_ref)
print('D-MnOOH: dH=',dh_D_MnOOH, 'eV  ',dh_D_MnOOH*kjmol, 'kj/mol')

dh_B_MnOH2=B_MnOH2-(2*o_ref2+Mn_metal2+2*h_ref)
print('B-MnOH2: dH=',dh_B_MnOH2, 'eV  ',dh_B_MnOH2*kjmol, 'kj/mol')

dh_D_MnOH2=D_MnOH2-(2*o_ref2+Mn_metal2+2*h_ref)
print('D-MnOH2: dH=',dh_D_MnOH2, 'eV  ',dh_D_MnOH2*kjmol, 'kj/mol   exp: ', dh_D_MnOH2_exp*kjmol, 'kjmol')

dh_Mn3O4=Mn3O4-(4*o_ref2+3*Mn_metal23)
print('Mn3O4: dH=',dh_Mn3O4, 'eV  ',dh_Mn3O4*kjmol, 'kj/mol   exp: ', dh_Mn3O4_exp*kjmol, 'kjmol')

dh_Mn2O3=Mn2O3-(3*o_ref2+2*Mn_metal3)
print('Mn2O3: dH=',dh_Mn2O3, 'eV  ',dh_Mn2O3*kjmol, 'kj/mol   exp: ', dh_Mn2O3_exp, 'kjmol')




################################
print('\nCalculate dGs now')

TdS_MnO=(S_MnO_exp-S_Mn_metal_exp-S_o2_gas_exp/2)*convert_dS
TdS_MnO2=(S_MnO2_solid_exp-S_Mn_metal_exp-S_o2_gas_exp)*convert_dS
dg_B_MnO2=dh_B_MnO2-TdS_MnO2
print('B-MnO2: dG= ' ,dg_B_MnO2, 'eV ', dg_B_MnO2*kjmol, 'kj/mol   exp: ', dg_B_MnO2_exp, 'kj/mol')
dg_R_MnO2=dh_R_MnO2-TdS_MnO2
print('R-MnO2: dG= ' ,dg_R_MnO2, 'eV ', dg_R_MnO2*kjmol, 'kj/mol   exp: ', dg_R_MnO2_exp, 'kj/mol')
dg_D_MnO2=dh_D_MnO2-TdS_MnO2
print('D-MnO2: dG= ' ,dg_D_MnO2, 'eV ', dg_D_MnO2*kjmol, 'kj/mol   exp: ', dg_D_MnO2_exp, 'kj/mol')
dg_E_MnO2=dh_E_MnO2-TdS_MnO2
print('E-MnO2: dG= ' ,dg_E_MnO2, 'eV ', dg_E_MnO2*kjmol, 'kj/mol')

 
TdS_MnOOH=(S_B_MnOOH_exp-S_Mn_metal_exp-S_o2_gas_exp-S_h2_gas_exp/2)*convert_dS
#TdS_MnOOH=(-S_Mn_metal_exp-S_o2_gas_exp-S_h2_gas_exp)*convert_dS
#TdS_MnOOH=(S_D_MnOH2_exp-S_Mn_metal_exp-S_o2_gas_exp-S_h2_gas_exp/2)*convert_dS
dg_B_MnOOH=dh_B_MnOOH-TdS_MnOOH
print('B-MnOOH: dG= ' ,dg_B_MnOOH, 'eV ', dg_B_MnOOH*kjmol, 'kj/mol   exp: ', dg_B_MnOOH_exp, 'kj/mol')
dg_D_MnOOH=dh_D_MnOOH-TdS_MnOOH
print('D-MnOOH: dG= ' ,dg_D_MnOOH, 'eV ', dg_D_MnOOH*kjmol, 'kj/mol   exp: ', dg_D_MnOOH_exp, 'kj/mol')

TdS_MnOH2=(S_D_MnOH2_exp-S_Mn_metal_exp-S_o2_gas_exp-S_h2_gas_exp)*convert_dS
dg_B_MnOH2=dh_B_MnOH2-TdS_MnOH2
print('B-MnOH2: dG= ' ,dg_B_MnOH2, 'eV ', dg_B_MnOH2*kjmol, 'kj/mol')
dg_D_MnOH2=dh_D_MnOH2-TdS_MnOH2
print('D-MnOH2: dG= ' ,dg_D_MnOH2, 'eV ', dg_D_MnOH2*kjmol, 'kj/mol   exp: ', dg_D_MnOH2_exp, 'kj/mol')

TdS_Mn3O4=(S_Mn3O4_exp-3*S_Mn_metal_exp-2*S_o2_gas_exp)*convert_dS
dg_Mn3O4=dh_Mn3O4-TdS_Mn3O4
print('Mn3O4: dG= ' ,dg_Mn3O4, 'eV ', dg_Mn3O4*kjmol, 'kj/mol  exp: ', dg_Mn3O4_exp, 'kjmol')

TdS_Mn2O3=(S_Mn2O3_exp-2*S_Mn_metal_exp-3*S_o2_gas_exp/2)*convert_dS
dg_Mn2O3=dh_Mn2O3-TdS_Mn2O3
print('Mn2O3: dG= ' ,dg_Mn2O3, 'eV ', dg_Mn2O3*kjmol, 'kj/mol  exp: ', dg_Mn2O3_exp, 'kjmol\n')

print('MnO: TdS= ', TdS_MnO, 'eV')
print('MnO2: TdS= ', TdS_MnO2, 'eV')
print('MnOOH: TdS= ', TdS_MnOOH, 'eV')
print('MnOH2: TdS= ', TdS_MnOH2, 'eV')
print('Mn3O4: TdS= ', TdS_Mn3O4, 'eV')
print('Mn2O3: TdS= ', TdS_Mn2O3, 'eV')


print('\n#############################the end ###################################')
