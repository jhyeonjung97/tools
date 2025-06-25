#!/usr/bin/env python
"""
NO3RR (Nitrate Reduction Reaction) 오버포텐셜 계산 및 시각화 스크립트
"""
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import pow
from pylab import *
import matplotlib
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.ticker import *
from scipy import *
import subprocess
import csv
from matplotlib import rc
import seaborn as sns
from matplotlib.colors import ListedColormap

fig_width_pt = 1.8*246.0
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size = [fig_width,fig_height]

font_size = 9 
tick_font_size = 8
xlabel_pad = 8
ylabel_pad = 18
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.labelsize'] = font_size
matplotlib.rcParams['legend.fontsize'] = font_size
matplotlib.rcParams['xtick.labelsize'] = tick_font_size
matplotlib.rcParams['ytick.labelsize'] = tick_font_size
matplotlib.rcParams['lines.linewidth'] = 1.

## gibbs free energy correction
ddgno2 = 0.0
ddgno2h = 0.0
dgnh2 = 0.0
dgnh3 = 0.0

no2_solv = 0.0
no2h_solv = 0.0
nh2_solv = 0.0
nh3_solv = 0.0

ddg_corr = [ddgno2, ddgno2h, dgnh2, dgnh3]
solv_corr = [no2_solv, no2h_solv, nh2_solv, nh3_solv]

no3rr_onset = 0.70

steps = ['NO2H*->NO2*', 'NH2*->NH3*']
## steps = ['NO3(g)->NO3*', 'NO3*->NO3H*', 'NO3H*->NO2H*', 'NO2H*->NO2*', 'NO2*->NH2*', 'NH2*->NH3*', 'NH3*->NH3(g)']

def overpotential_no3rr_full(dgno, dgnoh, dgnh2, dgnh3):
    dg = [dgnoh-dgno, dgnh3-dgnh2]
    m = max(dg)
    return [round(m+no3rr_onset,2), -round(m,2), no3rr_step(dg.index(m))]

def no3rr_step(i):
    return str(steps[i])

surfs = ['Ni(111)', 'Cu(100)', 'Cu(111)', 'Ru(0001)', 'Pt(111)', 'Rh(111)', 'Ag(111)']
comps = ['Ni27', 'Cu36', 'Cu36', 'Ru36', 'Pt27', 'Rh27', 'Ag27']
facets = ['111', '100', '111', '0001', '111', '111', '111']
markers = {'111': 'o', '100': 's', '0001': 'd'}
colors = {'Ni27': 'dodgerblue', 'Cu36': 'gold', 'Ru36': 'darkorange', 'Pt27': 'green', 'Rh27': 'orangered', 'Ag27': 'mediumorchid'}

calc_systems = []
cathub_nitrate = pd.read_csv('~/bin/tools/fwp/cathub-nitrate.csv')
for i in range(7):
    comp = comps[i]
    facet = facets[i]
    cathub_nitrate_comp = cathub_nitrate[(cathub_nitrate['chemicalComposition'] == comp) & (cathub_nitrate['facet'] == facet)]
    NH3_ads = cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "2.5H2(g) + NO(g) + * -> NH3* + H2O(g)")]['reactionEnergy'].values[0]
    NH2_ads = cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "2.0H2(g) + NO(g) + * -> NH2* + H2O(g)")]['reactionEnergy'].values[0]
    OH_ads = -cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "OH* + 0.5H2(g) -> H2O(g) + *")]['reactionEnergy'].values[0]
    O_ads = OH_ads -cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "0.5H2(g) + O* -> OH*")]['reactionEnergy'].values[0]
    NO_ads = cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "NO(g) + * -> NO*")]['reactionEnergy'].values[0]
    if comp == 'Ag27':
        NOH_ads = cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "0.5H2(g) + NO(g) + * -> NHO*")]['reactionEnergy'].values[0]
    else:
        NOH_ads = cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "0.5H2(g) + NO(g) + * -> NOH*")]['reactionEnergy'].values[0]
    # NO2_ads = NO_ads + O_ads - cathub_nitrate_comp[(cathub_nitrate_comp['Equation'] == "NO2* + * -> NO* + O*")]['reactionEnergy'].values[0]
    # protonation = NOH_ads - NO_ads
    # NO2H_ads = NO2_ads + protonation
    calc_systems.append([NO_ads, NOH_ads, NH2_ads, NH3_ads, 0.0, surfs[i], colors[comp], markers[facet]])

# CSV file
with open('Final_dGs_for_M_system_NO3RR.csv', 'w') as myfile:
    fieldnames = ['Surface name', 'dGNO', 'dGNOH', 'dGNH2', 'dGNH3', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()
    
    for i in range(len(calc_systems)):
        for k in range(len(steps)-1):
            calc_systems[i][k] += solv_corr[k] + ddg_corr[k]

        recalculated_over = overpotential_no3rr_full(calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3])
        calc_systems[i][4] = recalculated_over[0]
        wr.writerow([calc_systems[i][5], 
                    calc_systems[i][0], 
                    calc_systems[i][1], 
                    calc_systems[i][2], 
                    calc_systems[i][3], 
                    calc_systems[i][4],
                    recalculated_over[0], 
                    recalculated_over[1], 
                    recalculated_over[2]])

# Scaling plot
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1, x2 = -3, 0
y1, y2 = -6, 0
ax.axis([x1, x2, y1, y2])

ax.set_xlabel(r'$\Delta$G$_{\sf NO}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf NO_xH_y}$ (eV)')

xdata = []
y1data = []
y2data = []
y3data = []
y4data = []
y5data = []

for i in range(len(calc_systems)):
    xdata.append(calc_systems[i][0])
    y1data.append(calc_systems[i][1])
    y2data.append(calc_systems[i][2])
    y3data.append(calc_systems[i][3])

fit1 = polyfit(xdata, y1data, 1)
fit_fn1 = poly1d(fit1)

fit2 = polyfit(xdata, y2data, 1)
fit_fn2 = poly1d(fit2)

fit3 = polyfit(xdata, y3data, 1)
fit_fn3 = poly1d(fit3)

print(fit_fn1, fit_fn2, fit_fn3)

delta = 0.01
xx = np.arange(x1, x2, delta)
ax.plot(xx, xx, '--', lw=1, dashes=(3,1), c='grey')
ax.plot(xx, fit_fn1[1]*xx+fit_fn1[0], '--', lw=1, dashes=(3,1), c='orangered', label='NOH* scaling')
ax.plot(xx, fit_fn2[1]*xx+fit_fn2[0], '--', lw=1, dashes=(3,1), c='green', label='NH2* scaling')
ax.plot(xx, fit_fn3[1]*xx+fit_fn3[0], '--', lw=1, dashes=(3,1), c='blue', label='NH3* scaling')

for i in range(len(calc_systems)): 
    ax.plot(calc_systems[i][0], calc_systems[i][0], mec='black', mfc=calc_systems[i][6], mew=0.8, zorder=4, marker=calc_systems[i][7], linestyle='', 
            label=calc_systems[i][5]+' (BEEF-vdW)'+' : %.2f V' %(-(calc_systems[i][4]+no3rr_onset)))
    ax.plot(calc_systems[i][0], calc_systems[i][1], marker=calc_systems[i][7], ms=3, color='orangered')
    ax.plot(calc_systems[i][0], calc_systems[i][2], marker=calc_systems[i][7], ms=3, color='green')
    ax.plot(calc_systems[i][0], calc_systems[i][3], marker=calc_systems[i][7], ms=3, color='blue')

ax.legend(bbox_to_anchor=(1.01, 1.02), loc='lower right', borderaxespad=0.5, ncol=2, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)
eq1 = r'$\Delta$G$_{\sf NOH}$=%.2f$\Delta$G$_{\sf NO}$%.2f' %(fit_fn1[1],fit_fn1[0])
eq2 = r'$\Delta$G$_{\sf NO2}$=%.2f$\Delta$G$_{\sf NO}$%.2f' %(fit_fn2[1],fit_fn2[0])
eq3 = r'$\Delta$G$_{\sf NO2H}$=%.2f$\Delta$G$_{\sf NO}$%.2f' %(fit_fn3[1],fit_fn3[0])

ax.text(-1.4, -2.1, eq1, fontsize='small', color='orangered')
ax.text(-1.4, -2.4, eq2, fontsize='small', color='green')
ax.text(-1.4, -2.7, eq3, fontsize='small', color='blue')

fig.savefig('NO3RR_scaling.pdf', bbox_inches='tight')
plt.close()

def dgnoh_dgno_scaling(dgno):
    dgnoh = fit_fn1[1]*dgno + fit_fn1[0]
    return dgnoh

def dgnh2_dgno_scaling(dgno):
    dgnh2 = fit_fn2[1]*dgno + fit_fn2[0]
    return dgnh2

def dgnh3_dgno_scaling(dgno):
    dgnh3 = fit_fn3[1]*dgno + fit_fn3[0]
    return dgnh3

def overpotential_from_scaling(x, dgno):
    dgnh2_from_scaling = dgnh2_dgno_scaling(dgno)
    dgnh3_from_scaling = dgnh3_dgno_scaling(dgno)
    dg = [x, dgnh3_from_scaling-dgnh2_from_scaling] ##
    m = max(dg)
    return m+no3rr_onset

# 1D plot
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

x1 = -2.0
x2 = 1.0
y2 = -1.0
y1 = 1.0

ax.axis([x1, x2, y1, y2])
xx = np.arange(x1, x2, delta)

ax.set_xlabel(r'$\Delta$G$_{\sf NOH}-\Delta$G$_{\sf NO}$ (eV)')
ax.set_ylabel(r'U$_{\sf NO3RR}$ (V)')
ax.set_ylim(ax.get_ylim()[::-1])

'''
NOH* scaling: 0.842 x - 0.1704  
NO2* scaling: 0.4684 x + 0.2301  
NO2H* scaling: 0.3104 x + 0.05975  
NH2* scaling: 0.4484 x - 3.415  
NH3* scaling: 0.28 x - 4.464

x = NO*
xx = NOH* - NO*
xx = (0.842 - 1)x + (-0.1704)
x = (xx + 0.1704) / (0.842 - 1)

step 1: NOH* - NO* -> xx
step 2: NO2H* - NO2* -> (0.3104-0.4684)x + (0.05975-0.2301)
    -> (0.3104-0.4684)*(xx+0.1704)/(0.842-1) + (0.05975-0.2301)
    -> (0.3104-0.4684)/(0.842-1)*xx + (0.3104-0.4684)*0.1704/(0.842-1) + (0.05975-0.2301)
step 3: NH3* - NH2* -> (0.28-0.4484)x + (-4.464-(-3.415))
    -> (0.28-0.4484)*(xx+0.1704)/(0.842-1) + (-4.464-(-3.415))
    -> (0.28-0.4484)/(0.842-1)*xx + (0.28-0.4484)*0.1704/(0.842-1) + (-4.464-(-3.415))

step 1: NOH* - NO* 
step 2: NO2H* - NO2* 
step 3: NH3* - NH2* 

'''

a1 = (0.3104-0.4684)/(0.842-1)
a2 = (0.3104-0.4684)*0.1704/(0.842-1) + (0.05975-0.2301)
b1 = (0.28-0.4484)/(0.842-1)
b2 = (0.28-0.4484)*0.1704/(0.842-1) + (-4.464-(-3.415))

print(f"a1 = {a1:.2f}, a2 = {a2:.2f}")
print(f"b1 = {b1:.2f}, b2 = {b2:.2f}")

plot(xx, -np.maximum(xx, 1.07*xx-0.87), '--', color='black', lw=0.67, dashes=(3,1), zorder=2)
plot(xx, no3rr_onset+0.0*xx, '--', color='black', lw=0.67, dashes=(3,1), zorder=2)

for i in range(len(calc_systems)):
    ax.plot(calc_systems[i][1]-calc_systems[i][0], -(calc_systems[i][4]-no3rr_onset),
            mec='black', mfc=calc_systems[i][6], mew=0.8, zorder=4,
            marker=calc_systems[i][7], linestyle='',
            label=calc_systems[i][5]+' (BEEF-vdW)'+' : %.2f V' %(-(calc_systems[i][4]+no3rr_onset)))

ax.legend(bbox_to_anchor=(1.0, 1.4), loc='upper right', borderaxespad=0.5, ncol=2, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)

fig.savefig('NO3RR_1D.pdf', bbox_inches='tight')
plt.close()

# 2D plot
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
zoom = 0.3
d1 = 2*zoom
d2 = 4*zoom
xcenter = 0.3
ycenter = -0.1

x1 = xcenter-d1
x2 = xcenter+d1
y1 = ycenter-d2
y2 = ycenter+d2
ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf NO2H*}$ - $\Delta$G$_{\sf NO2*}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf NO2*}$ (eV)')

x = np.arange(x1,x2+delta, delta)
y = np.arange(y1,y2+delta, delta)
X, Y = np.meshgrid(x, y)

Z = []
for j in y:
    tmp = []
    for i in x:
        tmp.append((overpotential_from_scaling(i,j)))
    Z.append(tmp)
Z = np.array(Z)

levels = np.arange(0.5, 2.0, 0.1)
CS = plt.contourf(X, Y, Z, levels,
                 cmap='RdBu',
                 extend='max',
                 origin='lower')

cbar = plt.colorbar(CS, ticks=list(np.arange(0.5, 2.0, step=0.1)))
cbar.ax.set_ylabel(r'$\eta_{\sf NO3RR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

for i in range(len(calc_systems)):
    ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],
            mec='black', mfc=calc_systems[i][6], mew=0.8, zorder=4,
            marker=calc_systems[i][7], linestyle='',
            label=calc_systems[i][5]+' (BEEF-vdW)'+' : %.2f V' % (calc_systems[i][4]))

ax.legend(bbox_to_anchor=(1.2, 1.5), loc='upper right', borderaxespad=0.5, ncol=2, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)

fig.savefig('NO3RR_contour.pdf', bbox_inches='tight')
plt.close()