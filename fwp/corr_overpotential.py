#!/usr/bin/env python
"""
CORR (CO2 Reduction Reaction) 오버포텐셜 계산 및 시각화 스크립트
"""
import matplotlib
import numpy as np
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

# 그래프 설정
fig_width_pt = 1.8*246.0
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size = [fig_width,fig_height]
fig = plt.figure(figsize=fig_size,dpi=300)

# 폰트 및 스타일 설정
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

# 축 설정
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
zoom = 0.3
d1 = 3*zoom
d2 = 4*zoom
xcenter = 0.5
ycenter = -1.0

x1 = xcenter-d1
x2 = xcenter+d1
y1 = ycenter-d2
y2 = ycenter+d2
ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf COH*}$ - $\Delta$G$_{\sf CO*}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf CO*}$ (eV)')

# 메시 그리드 생성
delta = 0.01
x = np.arange(x1,x2+delta, delta)
y = np.arange(y1,y2+delta, delta)
X, Y = np.meshgrid(x, y)

# 스케일링 함수
def dgc_dgco_scaling(dgco):
    dgc = 2.2*dgco + 2.34 
    return dgc

# 상수 정의
corr_onset = 0.26
ch4 = -2.16

# 오버포텐셜 계산 함수
def overpotential_from_dgc_dgco_scaling(x, dgco):
    dgc_from_scaling = dgc_dgco_scaling(dgco)
    dg14 = [dgco, x, dgc_from_scaling -(x+dgco), (ch4-dgc_from_scaling)/4.0]
    m = max(dg14)
    return m+corr_onset

def overpotential_corr_full(dgco, dgcoh, dgc, dgch):
    dg23 = [dgco, dgcoh-dgco, dgc-dgcoh, dgch-dgc, (ch4-dgch)/3.0]
    m = max(dg23)
    return [round(m+corr_onset,2), -round(m,2), corr_step(dg23.index(m))]

# 단계 정의
steps = ['CO(g)+*->CO*','CO*->COH*', 'COH*->C*', 'C*->CH*', 'CH*->CH4(g)']
def corr_step(i):
    return str(steps[i])

# Z 값 계산
Z = []
for j in y:
    tmp = []
    for i in x:
        tmp.append((overpotential_from_dgc_dgco_scaling(i,j)))
    Z.append(tmp)
Z = np.array(Z)

# 컨투어 플롯 생성
levels = np.arange(0.0, 1.5, 0.1)
CS = plt.contourf(X, Y, Z, levels,
                 cmap=ListedColormap([
                    '#a50026', '#d73027', '#f46d43', '#fdae61',
                    '#fee090', '#ffffbf', '#ffffe5', '#e0f3f8',
                    '#abd9e9', '#74add1', '#4575b4', '#313695',
                    '#d0d1e6'
                 ]),
                 extend='max',
                 origin='lower')

# 컬러바 설정
cbar = plt.colorbar(CS, ticks=list(np.arange(0.0, 1.4, step=0.1)))
cbar.ax.set_ylabel(r'$\eta_{\sf CORR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# 축 설정
ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# 계산 시스템 데이터
calc_systems = []
colors = ['#a6cee3','#fdbf6f','#b2df8a','#1f78b4','#e31a1c','#fb9a99','#33a02c']
tcolor = '#d9d9d9'
bcolor = '#252525'
symbols = ['D','^','o','s','h','H','x','*','+']

# 시스템 데이터 추가
ddgco = 0.09+0.45
ddgcoh = 0.37+0.485
ddgc = 0.08
ddgch = 0.33

co_solv = -0.06
coh_solv = -0.1
coh_solv_bonus = -0.3
c_solv = 0.0
ch_solv = 0.0

hc_dest = 0.65

# 시스템 데이터 추가
calc_systems.append([-1.76+hc_dest, -1.47+hc_dest, -1.71, -1.94, 0.0, r'Rh(100) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'grey', 'o', 'purple'])
calc_systems.append([-1.75+hc_dest, -1.20+hc_dest, -1.26, -1.26, 0.0, r'Rh(211) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'grey', 'D', 'purple'])
calc_systems.append([-1.63+hc_dest, -0.94+hc_dest, -1.28, -1.27, 0.0, r'Pd(100) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'brown', 'o', 'brown'])
calc_systems.append([-1.65+hc_dest, -0.96+hc_dest, -1.16, -1.16, 0.0, r'Pd(211) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'brown', 'D', 'brown'])
calc_systems.append([-1.65+hc_dest, -1.08+hc_dest, -1.06, -1.54, 0.0, r'Pt(100) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'green', 'o', 'darkblue'])
calc_systems.append([-1.71+hc_dest, -1.14+hc_dest, -0.52, -0.52, 0.0, r'Pt(211) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'green', 'D', 'darkblue'])
calc_systems.append([-1.64+hc_dest, -1.34+hc_dest, -1.66, -2.24, 0.0, r'Ni(100) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'green', 'o', 'green'])
calc_systems.append([-1.70+hc_dest, -1.04+hc_dest, -1.30, -1.30, 0.0, r'Ni(210) (0.66 ML CO*)', 'black', 0.0, 0.08, 1.5, '', 'green', 'D', 'green'])
calc_systems.append([-0.65, 0.34+coh_solv_bonus, 0.38, -0.39, 0.0, r'Cu(100) (Beef, ++solv.)', 'black', 0.0, 0.08, 1.5, '', 'lime', 'o', 'red'])
calc_systems.append([-0.65, 0.3+0.4, 0.38, -0.39, 0.0, r'Cu(100) (Beef, ++solv.)', 'black', 0.0, 0.08, 1.5, '', 'lime', 'o', 'orange'])
calc_systems.append([-0.49, 0.62+coh_solv_bonus, 1.51, 0.14, 0.0, r'Cu(111) (++solv.)', 'black', 0.0, 0.08, 1.5, '', 'lime', '^', 'red'])
calc_systems.append([-0.6, 1.10+coh_solv_bonus, 1.51, 0.14, 0.0, r'Cu(111) (DMC, ++solv.)', 'black', 0.0, 0.08, 1.5, '', 'lime', 'd', 'orangered'])
calc_systems.append([-1.21+hc_dest, 0.86+coh_solv_bonus, 0.90, -0.04, 0.0, r'Cu(211) (0.66 ML CO ++solv.)', 'black', 0.0, 0.08, 1.5, '', 'lime', 'D', 'red'])
calc_systems.append([-0.05, 1.37+coh_solv_bonus, 2.35, 0.99, 0.0, r'Au(111) (++solv.)', 'black', 0.0, 0.08, 1.5, '', 'gold', '^', 'gold'])

# CSV 파일로 결과 저장
with open('Final_dGs_for_MOOHx_system_CORR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dGCO','dGCOH','dGC','dGCH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    solv_corr = [co_solv, coh_solv, c_solv, ch_solv]
    ddg_corr = [ddgco, ddgcoh, ddgc, ddgch]
    
    for i in range(len(calc_systems)):
        for k in range(len(steps)-1):
            calc_systems[i][k] += solv_corr[k] + ddg_corr[k]

        recalculated_over = overpotential_corr_full(calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3])
        calc_systems[i][4] = recalculated_over[0]
        wr.writerow([calc_systems[i][5], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3], 
                    recalculated_over[0], recalculated_over[1], recalculated_over[2]])

# dE 값으로 CSV 파일 생성
with open('Final_dEs_for_MOOHx_system_CORR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dECO','dECOH','dEC', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print('Final dEs for MOOHx_system system')

    for i in range(len(calc_systems)):
        recalculated_over = overpotential_corr_full(calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3])
        calc_systems[i][4] = recalculated_over[0]
        wr.writerow([calc_systems[i][5], 
                    calc_systems[i][0]-ddgco,
                    calc_systems[i][1]-ddgcoh,
                    calc_systems[i][2]-ddgc,
                    calc_systems[i][3]-ddgch,
                    recalculated_over[0],
                    recalculated_over[1],
                    recalculated_over[2]])

# 플롯 생성
for i in range(len(calc_systems)):
    ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],
            mec=calc_systems[i][6], mfc=calc_systems[i][13], mew=0.8, zorder=4,
            marker=calc_systems[i][12], linestyle='',
            label=calc_systems[i][5]+' : %.2f V' % (calc_systems[i][4]))

# 범례 추가
ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=1, ncol=3, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)

# 그래프 저장
fig.savefig('CORR_contour.pdf', bbox_inches='tight')
plt.close()

# 스케일링 플롯 생성
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1 = -1.4
x2 = 0.1
ax.axis([x1, x2, x1, dgc_dgco_scaling(x2)])

ax.set_xlabel(r'$\Delta$G$_{\sf CO}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf COH}$,$\Delta$G$_{\sf C}$ (eV)')

xdata = []
ydata = []
y2data = []

for i in range(len(calc_systems)):
    xdata.append(calc_systems[i][0])
    ydata.append(calc_systems[i][2])
    y2data.append(calc_systems[i][1])

# 다항식 피팅
fit = polyfit(xdata, ydata, 1)
fit_fn = poly1d(fit)
print(fit_fn)
aa = fit_fn[1]
bb = fit_fn[0]

fit1 = polyfit(xdata, y2data, 1)
fit_fn1 = poly1d(fit1)
print(fit_fn1)

# 스케일링 오차 계산
for i in range(len(calc_systems)):
    error = calc_systems[i][2]-(fit_fn[1]*calc_systems[i][0]+fit_fn[0])
    # print(error, calc_systems[i])

xx = np.arange(x1, x2, delta)
ax.plot(xx, fit_fn[1]*xx+fit_fn[0], '--', lw=1, dashes=(3,1), c='black', label='C* scaling')
ax.plot(xx, xx, '--', lw=1, dashes=(3,1), c='grey')
ax.plot(xx, fit_fn1[1]*xx+fit_fn1[0], '--', lw=1, dashes=(3,1), c='orangered', label='COH* scaling')

# 데이터 포인트 플롯
for i in range(len(calc_systems)): 
    ax.plot(calc_systems[i][0], calc_systems[i][2],
           marker=calc_systems[i][12], ms=3,
           color='black')
    ax.plot(calc_systems[i][0], calc_systems[i][1],
           marker=calc_systems[i][12], ms=3,
           color='orangered')
    ax.plot(calc_systems[i][0], calc_systems[i][3],
           marker=calc_systems[i][12], ms=3,
           color='gold')
    ax.plot(calc_systems[i][0], calc_systems[i][0],
           mec=calc_systems[i][6], mfc=calc_systems[i][13], mew=0.8, zorder=4,
           marker=calc_systems[i][12], linestyle='',
           label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]+corr_onset)))

# 범례 및 방정식 추가
ax.legend(bbox_to_anchor=(-0.15, 1.55), loc=2, borderaxespad=0.5, ncol=3, 
         fancybox=True, shadow=True, fontsize='x-small', handlelength=2)
eq1 = r'$\Delta$G$_{\sf C}$=%.2f$\Delta$G$_{\sf CO}$+%.2f' %(fit_fn[1],fit_fn[0])
ax.text(-1.1, 2.5, eq1, fontsize='small', color='black')
eq3 = r'$\Delta$G$_{\sf COH}$=%.2f$\Delta$G$_{\sf CO}$+%.2f' %(fit_fn1[1],fit_fn1[0])
ax.text(-1.1, 2.2, eq3, fontsize='small', color='orangered')
eq2 = r'$\Delta$G$_{\sf CH}$=xx$\Delta$G$_{\sf CO}$+xx'
ax.text(-1.1, 1.9, eq2, fontsize='small', color='gold')

fig.savefig('CORR_scaling.pdf', bbox_inches='tight')
plt.close()

# 1D 플롯 생성
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

x1 = -1.0
x2 = 2.0
y2 = -2.0
y1 = -0.0

ax.axis([x1, x2, y1, y2])
delta = 0.01
x = np.arange(x1, x2, delta)

ax.set_xlabel(r'$\Delta$G$_{\sf COH}-\Delta$G$_{\sf CO}$ (eV)')
ax.set_ylabel(r'U$_{\sf CORR}$ (V)')
ax.set_ylim(ax.get_ylim()[::-1])

co_fixed = -0.5
plot(x, -np.maximum(x, -1.5 +2.34 -x), '--', color='black', lw=0.67, dashes=(3,1), zorder=2)
plot(x, corr_onset+0.0*x, '--', color='black', lw=0.67, dashes=(3,1), zorder=2)

# 데이터 포인트 플롯
for i in range(len(calc_systems)):
    ax.plot(calc_systems[i][1]-calc_systems[i][0], -(calc_systems[i][4]-corr_onset),
            mec=calc_systems[i][6], mfc=calc_systems[i][13], mew=0.8, zorder=4,
            marker=calc_systems[i][12], linestyle='',
            label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]-corr_onset)))

# 범례 추가
ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=0.5, ncol=3, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)

fig.savefig('CORR_1D.pdf', bbox_inches='tight')
plt.close() 