from sys import exit
import math
import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel, RationalQuadratic
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_validate
from sklearn.pipeline import Pipeline
from sklearn.utils import shuffle
from sklearn.metrics import mean_absolute_error as mae
from sklearn.feature_selection import SequentialFeatureSelector as SFS
from sklearn.feature_selection import  SelectFromModel
from sklearn.model_selection import KFold
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import seaborn


from mendeleev import element

import ase
import pylab as p
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 10}

matplotlib.rc('font', **font)

""" Colors"""


c_bulk = {2:'#d62728',
          3: '#ff7f0e',
          4: '#2ca02c',
          5: '#1f77b4',
          6: '#9467bd'}

""" Define dataset"""

d3 = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
d4 = ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn']
d5 = ['La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb']

include_d0 = False
plot_d0 = False

outliers_remove = [#('W', 3),
                   #('Os', 3),
                   ('Sc', 2),
                   ('Cu', 2),
                   ('Mn', 6),
                   ('')
                   #('Ni', 2)
                   ]
all_data = pd.read_csv('create_data/all_data_dos_up_down.csv', index_col=0)



#all_data_old = pd.read_csv('all_data.csv', index_col=0)
#print(all_data['O 2p-center'])
#print(all_data['Bulk ICOHP'].values - all_data_old['Bulk ICOHP'].values)
#exit()
idx_remove = []
for i, dat in all_data.iterrows():
    if (dat['Metal Symbol'], dat['Oxidation State']) in outliers_remove:
        idx_remove += [i]

all_data.drop(idx_remove, inplace=True)

all_data = all_data.rename(columns={"Oxidation State": "Oxidation state bulk",
                            "Surface Oxidation State": "Oxidation state surface"})


all_data_O = all_data[all_data['Adsorbate'] == 'O']
all_data_OH = all_data[all_data['Adsorbate'] == 'OH']
all_data_O_OH = all_data[all_data['Adsorbate'] == 'O-OH']

all_data = all_data_O.merge(all_data_OH, on=['Bulk ICOHP', 'Oxidation state bulk',
       'Metal Symbol', 'Outer Electrons', 'Bulk File',
       'Psi', 'O 2p-center', 'predicted COHP', 'row', 'p ICOHP','Oxidation state surface'],
                                  suffixes=(' O', ' OH'))

all_data = all_data.merge(all_data_O_OH, on=['Bulk ICOHP', 'Oxidation state bulk',
       'Metal Symbol', 'Outer Electrons', 'Bulk File',
       'Psi', 'O 2p-center', 'predicted COHP', 'row', 'p ICOHP','Oxidation state surface'])

all_data = all_data.rename(columns={"Adsorption Energy": "Adsorption Energy O-OH"})

 # & (all_data['Oxidation state']==4)]

#def overpot(O_OH):
#    #O_OH += 0.001 - 0.339
#    return max([O_OH, 3.2-O_OH]) - 1.23

def get_overpotential(dGOH, dGO, dGOOH):
    """
    ideal step = 1.23 eV

                         _O2_
                   _OOH_|  ->  step4
               _O_|    ->  step3
          _OH_|   ->  step2
    _H2O_|  ->  step1
    """
    step_1 = dGOH
    step_2 = dGO - dGOH
    step_3 = dGOOH - dGO
    step_4 = 4 * 1.23 - dGOOH
    #step_energies = [step1, step2, step3, step4]
    step = []
    overpotential = []
    for i in range(len(dGOH)):


        #overpotential += [np.max([step_1[i], step_2[i], step_3[i], step_4[i]]) - 1.23]
        overpotential += [np.max([step_1[i], step_2[i], step_3[i], step_4[i]]) - 1.23]
    #step_energies = [step2, step3]  # Only consider step2 and step 3 for simplicity
    #overpotential = np.max(step_energies)
    return overpotential#, step

def ooh_oh_scaling(doh):
    """Oxide scaling for *OH vs *OOH (from previous studies)"""
    dooh = 0.8443 * doh + 3.136
    return dooh

def get_volcano(x):
    a = 1.6299  # O vs OH
    b = 0.9523  # O vs OH
    c = 0.8443 # OOH vs OH
    d = 3.136 # OOH vs OH

    step_1 = (x-b) / (a-1)  # H2O -> OH step\n",
    step_2 = x  # OH -> O step
    step_3 = (c-a) * (x-b) / (a-1)  + d - b  # O -> OOH step
    step_4 = 4 * 1.23 - c * (x+b) / (a-1) - d  # OOH -> O2 step\n",

    volcano_line = []
    for i in range(len(x)):
        print(np.argmax([step_1[i], step_2[i], step_3[i], step_4[i]]))
        volcano_line += [np.max([step_1[i], step_2[i], step_3[i], step_4[i]]) - 1.23]

    return volcano_line


outer_electrons_surface = []
outer_electrons_bulk = []
ionenergies_surface = []
ionenergies_bulk = []
is_d0 = []
group_nos = []
period_nos = []
atom_nos = []
psi = []
electronegativity = []
"""
if metal_symbol in ['Au', 'Ag']:
    #psi = get_outer_electrons(metal_symbol, oxidation_state) ** 2/element(metal_symbol).en_pauling ** 0.5
    #num = (int(element(metal_symbol).group.group_id) * 6**5) ** (1/3)
    #denom = (element(metal_symbol).en_pauling*element('O').en_pauling **5 ) ** (1/6)
    psi = num/denom
    #psi = int(element(metal_symbol).group.group_id) **c 2/element(metal_symbol).en_pauling ** 0.5
else:
    #psi = get_outer_electrons(metal_symbol, oxidation_state) ** 2/element(metal_symbol).en_pauling
    #psi = int(element(metal_symbol).group.group_id) ** 2/element(metal_symbol).en_pauling
    psi = num/denom

"""
for i, row in all_data.iterrows():
    m = row['Metal Symbol']
    ox_surf = float(row['Oxidation state surface'])
    ox_bulk = row['Oxidation state bulk']
    #for d in [d3, d4, d5]:
    #    if m in d:
    #m =

    Ele = element(m)
    group_no = int(Ele.group.group_id)
    period_nos += [Ele.period]
    atom_nos += [Ele.atomic_number]
    n_d = group_no - ox_surf
    n_d_b = group_no - ox_bulk
    group_nos += [group_no]
    #d_valence += [id + 3]
    outer_electrons_surface += [n_d]
    outer_electrons_bulk += [n_d_b]

    electronegativity += [Ele.en_pauling]
    #num = group_no #( n_d * 6 **4) ** (2/5)  #n_d**2   #(n_d * 6**4) ** (2/5)
    #denom = (Ele.en_pauling * element('O').en_pauling**5) **(1/6)   # * element('O').en_pauling**4) ** (1/5)
    #psi += [num/denom]

    if n_d_b == 0:
        is_d0 += [1]
    else:
        is_d0 += [0]
    ions = element(m).ionenergies
    c = ox_surf - math.floor(ox_surf)
    ox_check = math.floor(ox_surf) * (1-c) + math.ceil(ox_surf) * c
    #assert ox_check == ox_surf
    #for ox in [0,1,2,3,4]:
    #    print(ions.get(ox))

    ionenergies_bulk += [ions.get(ox_bulk+1) / ox_bulk]  # ox_bulk
    ionenergies_surface += [(ions.get(math.floor(ox_surf+1)) * (1-c) + ions.get(math.ceil(ox_surf+1)) * c)/ox_surf] #
    #ionenergies_surface += [ions.get(math.ceil(ox_surf+1)) / math.ceil(ox_surf)] #
    #icohp_surface += [(ions.get(math.floor(ox_surf)+1) * (1-c) + ions.get(math.ceil(ox_surf)+1) * c)/ox_surf]


all_data['Outer electrons surface'] = outer_electrons_surface
all_data['Outer electrons bulk'] = outer_electrons_bulk
all_data['Is d0'] = is_d0

all_data['Normalized ICOHP'] = all_data['Bulk ICOHP'] / all_data['Oxidation state bulk']
all_data['Normalized ICOHP Pred'] = all_data['predicted COHP'] / all_data['Oxidation state bulk']

all_data['Ionization energy surface'] = ionenergies_surface
all_data['Ionization energy bulk'] = ionenergies_bulk

all_data['Electronegativity'] = electronegativity

#all_data['Psi'] = psi

all_data['Group'] = group_nos
all_data['Period'] = period_nos
all_data['Atomic nos'] = atom_nos



d0_data = all_data[all_data['Is d0'] == 1]
if not include_d0:
    all_data = all_data[all_data['Is d0'] == 0]


with open('all_data_backup.csv', 'w') as f:
    all_data.to_csv(f)


""" Plot ICOHP pred vs ICOHP"""
plot_icohp_pred = False
if plot_icohp_pred:
    from sklearn.metrics import r2_score
    import matplotlib.cm as cm
    xaxis = 'Normalized ICOHP'
    yaxis = 'Normalized ICOHP Pred'
    p.figure(figsize=(4.5,4.5))
    x = all_data['Oxidation state bulk']
    colors_array = [c_bulk[ox] for ox in x]
    p.scatter(all_data[xaxis].values, all_data[yaxis].values, color=colors_array,  marker='o')

    for ox_i in [2,3,4,5,6]:
        c = c_bulk[ox_i]
        p.scatter(0,0, color=c, marker='o', label='+{}'.format(ox_i))

    p.xlabel('Calculated normalized ICOHP (eV)', fontsize=12)
    p.ylabel('Predicted normalized ICOHP (eV)', fontsize=12)

    MAE = mae(all_data[xaxis].values, all_data[yaxis].values)
    R2 = r2_score(all_data[xaxis].values, all_data[yaxis].values)
    p.plot(all_data[xaxis].values, all_data[xaxis].values, '--', color='gray', label='MAE={}eV, R$^2$={}'.format(round(MAE,3), round(R2,2)))

    p.legend()
    p.xlim(-4.2,-0.5)
    p.ylim(-4.2,-0.5)
    #p.yticks([-15, -10, -5, 0])
    border = 0.15
    p.subplots_adjust(top=0.95, bottom=border, left=border, right=0.95)
    p.savefig('../overleaf_copy/figures/Fig_SI_ICOHP_pred.pdf')
    p.show()
    #exit()


#import sys
#sys.exit()
#


adsorbate = 'OH'
#all_data_ads = all_data[all_data['Adsorbate'] == adsorbate]

for x_axis in []:#['O 2p-center', 'Normalized ICOHP', 'Ionization energy bulk']:
    p.figure(figsize=(4.5,4.5))

    y_axis = 'Adsorption Energy {}'.format(adsorbate)

    for i, ox in enumerate([2, 3, 4, 5, 6]):
        color = c_bulk[ox]

        dat = all_data[(all_data['Oxidation state bulk'] == ox)] # & (all_data['Oxidation state surface'] == 4)]

        #metal_symbols = dat['Metal Symbol'].values
        #id0 = dat['Is d0'].values
        #color_array = [color if not m in ['Ga','Ge', 'Sn', 'In', 'Pb', 'Tl'] and not id0[i] else 'gray' for i,m in enumerate(metal_symbols)]
        p.scatter(dat[x_axis], dat[y_axis], marker='o', color=color,#color_array,
                label = '+{}'.format(ox))
        #color_array,
                #label = '+{}'.format(ox))
        #for i in range(len(dat)):
        #    p.text(dat[x_axis].values[i], dat[y_axis].values[i], dat['Metal Symbol'].values[i])
    #datd0 = d0_data[(d0_data['Oxidation state bulk'] == ox)]
    if plot_d0:
        p.scatter(d0_data[x_axis], d0_data[y_axis], marker='x', color='gray', label='d$^0$')

    X = all_data[x_axis].values
    Y = all_data[y_axis].values
    #X = dat[x_axis].values
    #Y = dat[y_axis].values
    fit = 1
    if x_axis == 'Normalized ICOHP':
        #def fit_func(x, ro, b, c):
        #    return -1 * np.exp((x-ro)/b) + c
        #popt, pcov = curve_fit(fit_func, X, Y)
        coef = np.polyfit(X, Y, 2)
        fit = np.poly1d(coef)
        X_seq = np.linspace(-4.1,-0.8,100).reshape(-1, 1)

        #X_seq = np.linspace(-4.1,0,100).reshape(-1, 1)

        MAE = mae(Y, fit(X))#np.mean(np.abs(errors))
        #MAE = mae(Y, fit_func(X, *popt))#np.mean(np.abs(errors))
        p.plot(X_seq, fit(X_seq) , "k--", alpha=0.9, label='MAE: {} eV'.format(MAE.round(2)))
        #p.plot(X_seq, fit_func(X_seq, *popt), "k--", alpha=0.9, label='MAE: {} eV'.format(MAE.round(2)))

    elif x_axis == 'O 2p-center':
        coef = np.polyfit(X, Y, 1)
        fit = np.poly1d(coef)
        X_seq = np.linspace(-6.5,-1,100).reshape(-1, 1)
        #errors = Y - fit(X)
        MAE = mae(Y, fit(X)) #np.mean(np.abs(errors))
        p.plot(X_seq, fit(X_seq) , "k--", alpha=0.9, label='MAE: {} eV'.format(MAE.round(2)))

    elif x_axis == 'Ionization energy bulk':
        #def fit_func(x, ro, b, c):
        #    return -1 * np.exp((x-ro)/b) + c
        #popt, pcov = curve_fit(fit_func, X, Y)
        coef = np.polyfit(X, Y, 2)
        fit = np.poly1d(coef)
        X_seq = np.linspace(min(X), 20,100).reshape(-1, 1)

        #X_seq = np.linspace(-4.1,0,100).reshape(-1, 1)

        MAE = mae(Y, fit(X))#np.mean(np.abs(errors))
        #MAE = mae(Y, fit_func(X, *popt))#np.mean(np.abs(errors))

        p.plot(X_seq, fit(X_seq) , "k--", alpha=0.9, label='MAE: {} eV'.format(MAE.round(2)))
                #p.plot(X_seq, fit_func(X_seq, *popt), "k--", alpha=0.9, label='MAE: {} eV'.format(MAE.round(2)))

    if adsorbate == 'OH':
        p.ylabel('$\Delta\mathrm{E_{OH}}$ (eV)', fontsize=12)
    if adsorbate == 'O':
        p.ylabel('$\Delta \mathrm{E_{O}}$ (eV)', fontsize=12)
    if adsorbate == 'O-OH':
        p.ylabel('$\Delta \mathrm{E_{O}} - \Delta \mathrm{E_{OH}}$ (eV)', fontsize=12)
    p.xlabel(x_axis + ' (eV)', fontsize=12)


    #p.hlines(4.65, -4, 1, linestyles='dashed', color='k', alpha=0.5)
    p.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
    p.gca().yaxis.set_major_locator(MultipleLocator(1))

    if x_axis == 'Normalized ICOHP':
        p.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
        p.gca().xaxis.set_major_locator(MultipleLocator(0.5))
    elif x_axis == 'Ionization energy bulk':
        p.gca().xaxis.set_minor_locator(MultipleLocator(0.5))
        p.gca().xaxis.set_major_locator(MultipleLocator(2))
    else:
        p.gca().xaxis.set_minor_locator(MultipleLocator(0.2))
        p.gca().xaxis.set_major_locator(MultipleLocator(1))

    #p.xlim(-3.65, 0.2)
    p.legend(loc='lower right')
    border = 0.15
    p.subplots_adjust(top=0.95, bottom=border, left=border, right=0.95)

    if plot_d0:
        p.savefig('../overleaf_copy/figures/Fig_{}_vs_ads{}-with-d0.pdf'.format(x_axis, adsorbate))
    else:
        p.savefig('../overleaf_copy/figures/Fig_{}_vs_ads{}.pdf'.format(x_axis, adsorbate))
    p.show()
    #exit()
#

correlation_matrix = 0

if correlation_matrix:
    all_features = [
        'Adsorption Energy O',
        'Adsorption Energy OH',
        'Adsorption Energy O-OH',
        'O 2p-center',
        'Normalized ICOHP',
        'Ionization energy surface',
        'Ionization energy bulk',
        'Group',
        'Outer electrons surface',
        'Outer electrons bulk',
        #'Electronegativity',
        'Oxidation state surface',
        'Oxidation state bulk',

        ]

    p.figure(figsize=(9,9))
    ax = p.gca()
    X = all_data[all_features]
    from sklearn.preprocessing import StandardScaler
    #scaler = StandardScaler()
    #X = scaler.fit_transform(X)
    #from numpy import corrcoef
    corr = X.corr() #np.corrcoef(X.T) #X
    seaborn.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool),
                    cmap=seaborn.diverging_palette(220, 10, as_cmap=True),
                    square=True, ax=ax)

    p.subplots_adjust(top=0.95, bottom=0.25, left=0.25, right=0.95)
    p.savefig('../overleaf_copy/figures/SI_correlation_matrix')
    p.show()
    exit()

pca = 0
if pca:
    from sklearn import datasets, decomposition
    all_features = [
        'O 2p-center',
        'Normalized ICOHP',
        #'Ionization energy surface',
        'Ionization energy bulk',
        'Group',
        'Outer electrons surface',
        #'Outer electrons bulk',
        #'Oxidation state surface',
        'Oxidation state bulk',
        ]
    p.figure(figsize=(9,9))
    ax = p.gca()
    X = all_data[all_features]
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    pca = decomposition.PCA(n_components=6)
    pca.fit(X)
    reduced_X = pca.transform(X)
    corr = pca.components_.T

    seaborn.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool),
                    cmap=seaborn.diverging_palette(220, 10, as_cmap=True),
                        square=True, ax=ax)
    ax.set_yticklabels(all_features, rotation=0)
    #ax.set_xticklabels(['Comp. 1', 'Comp. 2', 'Comp. 3', 'Comp. 4'], rotation=90)
    p.subplots_adjust(top=0.95, bottom=0.25, left=0.25, right=0.95)
    #p.savefig('SI_correlation_matrix')
    p.show()
    #exit()

    print(pca.components_)
    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)
    #print(reduced_X)
    print(reduced_X.shape)

    y = all_data['Adsorption Energy OH'.format(adsorbate)].values
    kernel = WhiteKernel() + RBF(length_scale=1) + RBF(length_scale=3)
    model = GPR(kernel) #alpha=0.03)#model)#, n_restarts_optimizer=5)#, alpha = opt_alpha)

    kf = KFold(n_splits=int(len(y)))

    j = 0
    k = 0
    figures = []
    Test_diff_list = []
    all_predictions = []
    all_u = []
    all_metals = []
    p.figure(figsize=(4.5,4.5))
    maes_test = []
    maes_train = []
    for train_index, test_index in kf.split(reduced_X):
        #print(test_index)
        #print(j)
        j+=1
        X_train, X_test = reduced_X[train_index], reduced_X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        pipe = Pipeline([('scaler', StandardScaler()),
                        ('model', model),
                                    ])
        pipe.fit(X_train, y_train)
        y1 = pipe.predict(X_test)
        y2 = pipe.predict(X_train)
        maes_test += [mae(y_test,y1)]
        maes_train += [mae(y_train,y2)]

    print('MAE test', np.mean(maes_test))
    print('MAE train', np.mean(maes_train))
    exit()
scaling =0

if scaling:
    dOH = all_data['Adsorption Energy OH'].values + 0.2975
    dO = all_data['Adsorption Energy O'].values - 0.014
    print(len(dOH), len(dO) )
    color_array = [c_bulk[o] for o in all_data['Oxidation state bulk'].values]
    p.figure(figsize=(10,10))
    p.scatter(dOH, dO, c=color_array)
    for i in range(len(dOH)):
        p.annotate(all_data['Metal Symbol'].values[i] + ' + ' + str(all_data['Oxidation state bulk'].values[i]), (dOH[i], dO[i]))

    #p.title('O-OH Scaling relation', fontsize=20)
    p.xlabel('$\Delta$G$_{OH}$ (eV)', fontsize=16)
    p.ylabel('$\Delta$G$_{O}$ (eV)', fontsize=16)


    fit = np.polyfit(dOH,
                     dO,
                     1)

    a = fit[0]
    b = fit[1]
    xfit = np.array([-3, 3])
    yfit = a * xfit + b
    yfit_r = a * dOH + b
    from sklearn.metrics import r2_score

    r2 = r2_score(dO, yfit_r)

    p.plot(xfit, yfit, '--', label='y={:.4f}x + {:.4f}.  R2={:.2f}'.format(a,b,r2))
    p.legend(prop={'size':18})
    p.savefig('../overleaf_copy/figures/O_OH_scaling.pdf')
    p.show()
    exit()

volcano = 0


if volcano:
    dOH = all_data['Adsorption Energy OH'].values + 0.2975
    dO = all_data['Adsorption Energy O'].values - 0.014
    x = dO - dOH

    color_array = [c_bulk[o] for o in all_data['Oxidation state bulk'].values]
    dOOH = ooh_oh_scaling(dOH)
    y = get_overpotential(dOH, dO, dOOH)
    #print(step)
    x_grid = np.arange(min(x), max(x), 0.01)
    v = get_volcano(x_grid)

    p.plot(x_grid, v, 'k--')
    p.scatter(x, y, color=color_array, marker='o', zorder=2)


    for i in range(len(x)):
        if y[i] < 0.6:
            #p.annotate(all_metals[i], (xv[i]+xshift, yv[i]))
            p.annotate(all_data['Metal Symbol'].values[i] + ' + ' + str(all_data['Oxidation state bulk'].values[i]), (x[i], y[i]))

    p.xlabel(r'$\Delta G_{O} - \Delta G_{OH}$ (eV)', fontsize=12)
    p.ylabel('Overpotential (V)', fontsize=12)
    p.legend()
    p.ylim(1.4,0.1)
    p.xlim(0.4,2.7)
    p.tight_layout()

    p.savefig('../overleaf_copy/figures/volcano_real.pdf')
    #p.show()

    #exit()

    #dGOH =

    #yv_real =  get_overpotential(dGOH, dGO, dGOOH)



feature_list = [
    'Normalized ICOHP',
#    ]
#"""
    'Ionization energy bulk',
    'Group',
    'Oxidation state bulk',
    'Outer electrons surface',
    'O 2p-center',
    'Outer electrons bulk',
    'Oxidation state surface',
    'Ionization energy surface',
]
#"""


feature_list_OH = [
    'Normalized ICOHP',
    'Ionization energy bulk',
    'Group',
    'Oxidation state bulk',
    'Outer electrons surface',
    'O 2p-center',
    ]

feature_select_OH = [
    'Normalized ICOHP',
    'Ionization energy bulk',
    'Group',
    'Oxidation state bulk',
    'Outer electrons surface',
    'O 2p-center',
    'Outer electrons bulk',
    'Oxidation state surface',
    'Ionization energy surface',
    ]

feature_list_O = [
    'Normalized ICOHP',
    'Outer electrons bulk',
    'Oxidation state bulk',
    'Ionization energy bulk',
    'Oxidation state surface',
    'O 2p-center',
    ]

feature_select_O = [
    'Normalized ICOHP',
    'Outer electrons bulk',
    'Oxidation state bulk',
    'Ionization energy bulk',
    'Oxidation state surface',
    'O 2p-center',
    'Ionization energy surface',
    'Group',
    'Outer electrons surface',
    ]

feature_list_O_OH = [
    'Normalized ICOHP',
    'Group',
    'Oxidation state bulk',
    ]

feature_select_O_OH = [
    'Normalized ICOHP',
    'Group',
    'Oxidation state bulk',
    'Outer electrons bulk',
    'O 2p-center',
    'Ionization energy surface',
    'Oxidation state surface',
    'Outer electrons surface',
    'Ionization energy bulk',
    ]


kernel = WhiteKernel() + RBF(length_scale=1) + RBF(length_scale=3)

feature_select = True#False
plot_only = True
selected_features = ''

for adsorbate in ['O-OH']:#['O', 'O-OH']:#, 'O-OH']:#, 'OH', 'O-OH']:
    #all_data_ads = all_data[all_data['Adsorbate'] == adsorbate]
    if not feature_select:
        if adsorbate == 'O':
            feature_list = feature_list_O
        elif adsorbate == 'OH':
            feature_list = feature_list_OH
        elif adsorbate == 'O-OH':
            feature_list = feature_list_O_OH
        feature_numbers = [len(feature_list)]
        selected_features = feature_list
    else:
        if plot_only:
            if adsorbate == 'O':
                feature_list = feature_select_O
            elif adsorbate == 'OH':
                feature_list = feature_select_OH
            elif adsorbate == 'O-OH':
                feature_list = feature_select_O_OH
        feature_numbers = [f + 1 for f in range(len(feature_list))]
        selected_features = []
        selected_features_idx = []
        selected_features_all = []
    MAE_test = []
    std_test = []
    MAE_train = []
    M_unc_test = []
    Test_diff_list_all = []


    descriptor = all_data[feature_list].values

    all_ads = all_data['Adsorption Energy {}'.format(adsorbate)].values
    all_metals = all_data['Metal Symbol'].values
    all_ox = all_data['Oxidation state bulk'].values
    X = descriptor
    y = all_ads

    #X, y, al_met, al_ox = shuffle(X, y, all_metals, all_ox)
    al_met = all_metals
    al_ox = all_ox

    if feature_select:
        N_f = len(feature_list)
        for n_features in feature_numbers:
            #if n_features == N_f:
                #for n in range(N_f):
                #    if not n in selected_features_idx:
                #        selected_features_idx += [n]
            #    selected_features = feature_list #np.array(feature_list)[selected_features_idx]
                #continue

            kf = KFold(n_splits=int(len(y)))
            #kf = KFold(n_splits=5)  # fast test

            # from itertools import combinations
            # N_f = X.shape[1]
            # print(n_features)
            # comb = combinations(range(N_f), n_features)
            candidate_features = []
            if plot_only:
                #for n in range(n_features):
                candidate_features += [list(range(n_features))]

            else:
                for n in range(N_f):
                    if not n in selected_features_idx:
                        candidate_features += [list(selected_features_idx) + [n]]

            print('candidates', candidate_features)

            mean_MAE_test = []
            mean_MAE_train = []
            for c in candidate_features:
                print('checking', np.array(feature_list)[c])
                X_select = X[:, c]
                Test_MAE_list_0 = []
                Train_MAE_list_0 = []
                for train_index, test_index in kf.split(X_select):
                    X_train, X_test = X_select[train_index], X_select[test_index]
                    y_train, y_test = y[train_index], y[test_index]
                    model = GPR(kernel)
                    pipe = Pipeline([('scaler', StandardScaler()),
                                      ('model', model),
                                            ])
                    pipe.fit(X_train, y_train)
                    tr_pred, tr_pred_std = pipe.predict(X_train, return_std=True)
                    te_pred, te_pred_std = pipe.predict(X_test, return_std=True)
                    Train_MAE_list_0 += [mae(tr_pred, y_train)]
                    Test_MAE_list_0 += [mae(te_pred, y_test)]
                print('   MAE:', np.mean(Test_MAE_list_0))
                mean_MAE_test += [np.mean(Test_MAE_list_0)]
                mean_MAE_train += [np.mean(Train_MAE_list_0)]

            MAE_test += [np.min(mean_MAE_test)]
            MAE_train += [np.min(mean_MAE_train)]

            idx = np.argmin(mean_MAE_test)
            #print(idx)
            selected_features_idx = candidate_features[idx]

            #print(features_idx)
            selected_features = np.array(feature_list)[selected_features_idx]
            print('Selected features N={} :'.format(n_features), selected_features, np.min(mean_MAE_test))
            #selected_features_all += [selected_features]
            #pipe = Pipeline([('scaler', StandardScaler()),
            #                ('feature_selection', SFS(model, n_features_to_select=n_features, direction="forward", cv=20)),
            #                        ])
            #pipe.fit(X, y)
            #selected_features = np.array(feature_list)[pipe['feature_selection'].get_support()]
            #X = pipe['feature_selection'].transform(X)

            #print(selected_features)

            #continue

            #exit()
        feature_list = selected_features

        p.figure(figsize=(4.5, 5))

        #p.violinplot(np.array(Test_diff_list_all).T, showextrema=False)

        p.plot(feature_numbers, MAE_test, 'o--', label='test')
        #p.fill_between(feature_numbers,np.array(MAE_test) - 0.5 * np.array(std_test), np.array(MAE_test) + 0.5 * np.array(std_test), label='test uncertainty', color='grey', alpha=0.2)
        p.plot(feature_numbers, MAE_train, 'o--', label='train')
        p.xlabel('Number of descriptors', fontsize=12)
        p.ylabel('{} model MAE(eV)'.format(adsorbate), fontsize=12)
        border=0.15
        p.legend()
        p.gca().yaxis.set_minor_locator(MultipleLocator(0.01))
        p.gca().yaxis.set_major_locator(MultipleLocator(0.05))
        p.subplots_adjust(top=0.6, bottom=0.15, left=0.15,right=0.95)
        p.gca().set_xticks(feature_numbers)
        p.xlim(0.5, 9+0.5)
        par1 = p.gca().twiny()

        par1.set_xticks(feature_numbers)
        par1.set_xticklabels(selected_features, rotation=90)
        par1.set_xlim(0.5, 9+0.5)
        p.savefig('../overleaf_copy/figures/N_descr_{}.pdf'.format(adsorbate))

        #p.show()
        #exit()
    else:
        descriptor = all_data[feature_list].values
        X = descriptor
        Test_MAE_list = []
        Train_MAE_list = []

        print('Number of data points:', len(y))

        kf = KFold(n_splits=int(len(y)))

        j = 0
        k = 0
        figures = []
        Test_diff_list = []
        all_predictions = []
        all_u = []
        all_metals = []
        p.figure(figsize=(4.5,4.5))


        for train_index, test_index in kf.split(X):
            #print(test_index)
            #print(j)
            j+=1
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            m_train, m_test = al_met[train_index], al_met[test_index]
            ox_train, ox_test = al_ox[train_index], al_ox[test_index]

            model = GPR(kernel) #alpha=0.03)#model)#, n_restarts_optimizer=5)#, alpha = opt_alpha)
            pipe = Pipeline([('scaler', StandardScaler()),
                            ('model', model),
                                        ])
            pipe.fit(X_train, y_train)

            #print(pipe['model'].kernel_)

            tr_pred, tr_pred_std = pipe.predict(X_train, return_std=True)
            te_pred, te_pred_std = pipe.predict(X_test, return_std=True)
            all_metals += [m_test[0] +  ' +' + str(ox_test[0])]
            all_predictions += [te_pred[0]]
            all_u += [te_pred_std[0]]
            Test_diff_list += [te_pred[0] - y_test[0]]
            #print(Test_diff_list)
            Train_MAE_list += [mae(tr_pred, y_train)]
            Test_MAE_list += [mae(te_pred, y_test)]

            #for i in range(len(y_test)):
            #    p.annotate(m_test[i] + ' +' + str(ox_test[i]), (y_test[i], te_pred[i]))
        color_array = [c_bulk[o] for o in al_ox]
        p.text(np.min(y + 0.4),  np.max(y), 'Features:' + '\n   -{}'.format('\n   -'.join(selected_features)),
               va='top')
        p.errorbar(y, all_predictions, all_u, marker='None', color='gray', linestyle='None', zorder=1)#, alpha=0.5)
        p.scatter(y, all_predictions, color=color_array, marker='o', zorder=2)
        if adsorbate == 'OH':
            p.xlabel('Calculated $\Delta\mathrm{E_{OH}}$ (eV)', fontsize=12)
            p.ylabel('Predicted $\Delta\mathrm{E_{OH}}$ (eV)', fontsize=12)
        elif adsorbate == 'O':
            p.xlabel('Calculated $\Delta\mathrm{E_{O}}$ (eV)', fontsize=12)
            p.ylabel('Predicted $\Delta\mathrm{E_{O}}$ (eV)', fontsize=12)
        else:
            p.xlabel('Calculated $\Delta\mathrm{E_{O}} - \Delta\mathrm{E_{OH}}$ (eV)', fontsize=12)
            p.ylabel('Predicted $\Delta\mathrm{E_{O}} - \Delta\mathrm{E_{OH}}$ (eV)', fontsize=12)
        p.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        p.gca().xaxis.set_minor_locator(MultipleLocator(0.2))
        p.plot([min(y), max(y)], [min(y), max(y)], 'k--')
        #p.annotate('\n'.join(), xy=(-2,6))
        #p.title(adsorbate)
        a = p.axes([.75, .2, .15, .25])#, axisbg='y')
        p.violinplot(Test_diff_list, showextrema=False)
        p.ylim(-0.75,0.75)
        p.xticks([])
        print('RMSE', round(np.mean(np.array(Test_MAE_list)**2)**0.5, 3))
        print('MAE', round(np.mean(Test_MAE_list), 3))

        p.title('MAE: {} eV'.format(round(np.mean(Test_MAE_list), 3)), fontsize=10)
        border = 0.15
        p.subplots_adjust(top=0.95, bottom=border, left=border, right=0.95)
        MAE_test += [np.mean(Test_MAE_list)]
        std_test += [np.std(Test_MAE_list)]
        MAE_train += [np.mean(Train_MAE_list)]
        M_unc_test += [np.mean(all_u)]
        Test_diff_list_all += [np.abs(Test_diff_list)]
        #print(Test_diff_list_all)
        p.ylabel('model error (eV)')
        #p.yticks([])
        p.savefig('../overleaf_copy/figures/Fig3_b_GP_{}_{}_features.pdf'.format(adsorbate, len(feature_list)))
        p.show()
        if adsorbate == 'O-OH':
            p.figure(figsize=(4.5, 4))
            plt_rng = np.arange(-0.5,2.7,0.01)
            volcano_line = get_volcano(plt_rng)
            p.plot(plt_rng, volcano_line, 'k--')
            #p.plot(plt_rng,[overpot(a) for a in plt_rng], 'r--')
            xv_real = np.array(y) - 0.3115

            #dGOH =

            #ads_OH = all_data['Adsorption Energy OH']].values
            #ads_O = all_data['Adsorption Energy O']].values
            dOH = all_data['Adsorption Energy OH'].values + 0.2975
            dO = all_data['Adsorption Energy O'].values - 0.014
            dOOH = ooh_oh_scaling(dOH)

            yv_real =  get_overpotential(dOH, dO, dOOH)
            #yv_real = get_volcano(np.array(y) - 0.3115)#

            xv = np.array(all_predictions) - 0.3115
            yv = get_volcano(np.array(all_predictions) - 0.3115)

            MAE = mae(yv, yv_real)
            V_error_list = np.array(yv_real) - np.array(yv)
            print('MAE volcano V', MAE )
            #xv = np.array(y) - 0.3115
            #yv = [-1 * overpot(a) for a in y]
            #print(xv.shape, yv.shape)
            p.scatter(xv, yv, color=color_array, marker='o', zorder=2)

            lower_v =  np.array(get_volcano(np.array(all_predictions) - 0.5* np.array(all_u))) - 0.3115
            upper_v =  np.array(get_volcano(np.array(all_predictions) + 0.5* np.array(all_u))) - 0.3115
            yerr = upper_v - lower_v
            for i in range(len(xv)):
                if yv[i] < 0.8:#lower_v[i] < 0.5:
                #if 'Ir' in all_metals[i] or 'Ru' in all_metals[i] or 'Rh' in all_metals[i]:
                    #if xv[i] < 1.6:
                    #    xshift = -0.1
                    #else:
                    #    xshift = 0.1
                    xshift=0
                    p.annotate(all_metals[i], (xv[i]+xshift, yv[i]))
            #ax9.errorbar(tr_pred, [overpot(a) for a in tr_pred], xerr=tr_pred_std,


            p.errorbar(xv, yv, #xerr=all_u, #Test_diff_list,
                       xerr=all_u,
                       #yerr=yerr,
                       lw=0,
                       color='k',
                       marker='o',
                       elinewidth=2,capsize=2,capthick=1,markersize=0, alpha=0.1)
            #ax9.errorbar(te_pred, [overpot(a) for a in te_pred], xerr=te_pred_std,
            #p.errorbar(te_pred, [-1 * overpot(a) for a in y_test], xerr=te_pred_std,
            #        lw=0, marker='o',elinewidth=2,capsize=2,capthick=1,markersize=5, label='test set')
            p.xlabel(r'Predicted $\Delta G_{O} - \Delta G_{OH}$ (eV)', fontsize=12)
            p.ylabel('Predicted overpotential (V)', fontsize=12)
            p.legend()
            p.ylim(1.4,0.1)
            p.xlim(0.4,2.7)
            p.tight_layout()

            #a = p.axes([.75, 0.75, .15, .2])#, axisbg='y')
            #p.violinplot(V_error_list, showextrema=True)
            #p.ylim(-0.75,0.75)
            #p.xticks([])
            #p.title('MAE: {} V'.format(round(MAE, 3)), fontsize=10)
            #p.ylabel('model error (V)')
            p.savefig('../overleaf_copy/figures/volcano_pred.pdf')
            p.show()
