#!/usr/bin/env python
"""
NO3RR (Nitrate Reduction Reaction) 오버포텐셜 계산 및 시각화 스크립트
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
ax.set_xlabel(r'$\Delta$G$_{\sf NO2*}$ - $\Delta$G$_{\sf NO3*}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf NO3*}$ (eV)')

# 메시 그리드 생성
delta = 0.01
x = np.arange(x1,x2+delta, delta)
y = np.arange(y1,y2+delta, delta)
X, Y = np.meshgrid(x, y)

# 스케일링 함수
def dgno2_dgno3_scaling(dgno3):
    dgno2 = 1.2*dgno3 + 0.34 
    return dgno2

# 상수 정의
no3rr_onset = 0.26
nh3 = -2.16

# 오버포텐셜 계산 함수
def overpotential_from_dgno2_dgno3_scaling(x, dgno3):
    dgno2_from_scaling = dgno2_dgno3_scaling(dgno3)
    dg14 = [dgno3, x, dgno2_from_scaling -(x+dgno3), (nh3-dgno2_from_scaling)/4.0]
    m = max(dg14)
    return m+no3rr_onset

def overpotential_no3rr_full(dgno3, dgno2, dgno, dgnh):
    dg23 = [dgno3, dgno2-dgno3, dgno-dgno2, dgnh-dgno, (nh3-dgnh)/3.0]
    m = max(dg23)
    return [round(m+no3rr_onset,2), -round(m,2), no3rr_step(dg23.index(m))]

# 단계 정의
steps = ['NO3(g)+*->NO3*','NO3*->NO2*+O*', 'NO2*->NO*+O*', 'NO*->N*+O*', 'N*->NH3(g)']
def no3rr_step(i):
    return str(steps[i])

# Z 값 계산
Z = []
for j in y:
    tmp = []
    for i in x:
        tmp.append((overpotential_from_dgno2_dgno3_scaling(i,j)))
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
cbar.ax.set_ylabel(r'$\eta_{\sf NO3RR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# 축 설정
ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# 계산 시스템 데이터
calc_systems = []

# Rh(111) 시스템
# NO3* + * -> NO2* + O*: -1.2395878427487332
# NO2* + * -> NO* + O*: -1.0895914773573168
# NO* + * -> N* + O*: -0.269270736942417
# H2(g) + NO(g) + * -> H2O(g) + N*: -3.3389313691295683
calc_systems.append([-1.79, -0.35, -0.27, -0.27, 0.0, r'Rh(111)', 'black', 0.0, 0.08, 1.5, '', 'grey', 'o', 'purple'])

# Cu(111) 시스템
# NO3* + * -> NO2* + O*: -0.6919226635291125
# NO2* + * -> NO* + O*: 0.3461081371860928
# NO* + * -> N* + O*: 0.5986936011322541
# H2(g) + NO(g) + * -> H2O(g) + N*: -1.7287077950604726
calc_systems.append([-1.52, -0.69, -0.60, -0.60, 0.0, r'Cu(111)', 'black', 0.0, 0.08, 1.5, '', 'lime', 'D', 'red'])

# Ag(111) 시스템
# NO3* + * -> NO2* + O*: 0.5252383719998761
# NO2* + * -> NO* + O*: 1.1825559089920716
# NO* + * -> N* + O*: 2.584226393941208
# H2(g) + NO(g) + * -> H2O(g) + N*: -0.3205507281527389
calc_systems.append([-0.98, 0.53, 0.53, 0.53, 0.0, r'Ag(111)', 'black', 0.0, 0.08, 1.5, '', 'gold', '^', 'gold'])

# Ru(0001) 시스템
# NO3* + * -> NO2* + O*: -1.9787671580997994
# NO2* + * -> NO* + O*: -1.7664259040393517
# NO* + * -> N* + O*: -0.8328600867389468
# H2(g) + NO(g) + * -> H2O(g) + N*: -3.5951534986561455
calc_systems.append([-1.77, -0.98, -0.49, -0.49, 0.0, r'Ru(0001)', 'black', 0.0, 0.08, 1.5, '', 'green', 's', 'green'])

# Ni(111) 시스템
# NO3* + * -> NO2* + O*: -1.6447447535028914
# NO2* + * -> NO* + O*: -1.8649398255511187
# NO* + * -> N* + O*: -0.48547410291212145
# H2(g) + NO(g) + * -> H2O(g) + N*: -3.4728181184036657
calc_systems.append([-1.64, -0.70, -0.49, -0.49, 0.0, r'Ni(111)', 'black', 0.0, 0.08, 1.5, '', 'brown', 'h', 'brown'])

# Pt(111) 시스템
# NO3* + * -> NO2* + O*: -0.02006962098676013
# NO2* + * -> NO* + O*: -0.3463837204835727
# NO* + * -> N* + O*: 0.522064027874876
# H2(g) + NO(g) + * -> H2O(g) + N*: -3.012794619844499
calc_systems.append([-1.04, -0.02, 0.52, 0.52, 0.0, r'Pt(111)', 'black', 0.0, 0.08, 1.5, '', 'blue', 'H', 'darkblue'])

# CSV 파일로 결과 저장
with open('Final_dGs_for_MOOHx_system_NO3RR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dGNO3','dGNO2','dGNO','dGNH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(len(calc_systems)):
        recalculated_over = overpotential_no3rr_full(calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3])
        calc_systems[i][4] = recalculated_over[0]
        wr.writerow([calc_systems[i][5], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3], 
                    recalculated_over[0], recalculated_over[1], recalculated_over[2]])

# dE 값으로 CSV 파일 생성
with open('Final_dEs_for_MOOHx_system_NO3RR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dENO3','dENO2','dENO','dENH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(len(calc_systems)):
        recalculated_over = overpotential_no3rr_full(calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], calc_systems[i][3])
        calc_systems[i][4] = recalculated_over[0]
        wr.writerow([calc_systems[i][5], 
                    calc_systems[i][0],
                    calc_systems[i][1],
                    calc_systems[i][2],
                    calc_systems[i][3],
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
fig.savefig('NO3RR_contour.pdf', bbox_inches='tight')
plt.close()

# 스케일링 플롯 생성
fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1 = -1.4
x2 = 0.1
ax.axis([x1, x2, x1, dgno2_dgno3_scaling(x2)])

ax.set_xlabel(r'$\Delta$G$_{\sf NO3}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf NO2}$,$\Delta$G$_{\sf NO}$ (eV)')

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
ax.plot(xx, fit_fn[1]*xx+fit_fn[0], '--', lw=1, dashes=(3,1), c='black', label='NO* scaling')
ax.plot(xx, xx, '--', lw=1, dashes=(3,1), c='grey')
ax.plot(xx, fit_fn1[1]*xx+fit_fn1[0], '--', lw=1, dashes=(3,1), c='orangered', label='NO2* scaling')

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
           label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]+no3rr_onset)))

# 범례 및 방정식 추가
ax.legend(bbox_to_anchor=(-0.15, 1.55), loc=2, borderaxespad=0.5, ncol=3, 
         fancybox=True, shadow=True, fontsize='x-small', handlelength=2)
eq1 = r'$\Delta$G$_{\sf NO}$=%.2f$\Delta$G$_{\sf NO3}$+%.2f' %(fit_fn[1],fit_fn[0])
ax.text(-1.1, 2.5, eq1, fontsize='small', color='black')
eq3 = r'$\Delta$G$_{\sf NO2}$=%.2f$\Delta$G$_{\sf NO3}$+%.2f' %(fit_fn1[1],fit_fn1[0])
ax.text(-1.1, 2.2, eq3, fontsize='small', color='orangered')
eq2 = r'$\Delta$G$_{\sf NH}$=xx$\Delta$G$_{\sf NO3}$+xx'
ax.text(-1.1, 1.9, eq2, fontsize='small', color='gold')

fig.savefig('NO3RR_scaling.pdf', bbox_inches='tight')
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

ax.set_xlabel(r'$\Delta$G$_{\sf NO2}-\Delta$G$_{\sf NO3}$ (eV)')
ax.set_ylabel(r'U$_{\sf NO3RR}$ (V)')
ax.set_ylim(ax.get_ylim()[::-1])

no3_fixed = -0.5
plot(x, -np.maximum(x, -1.5 +0.34 -x), '--', color='black', lw=0.67, dashes=(3,1), zorder=2)
plot(x, no3rr_onset+0.0*x, '--', color='black', lw=0.67, dashes=(3,1), zorder=2)

# 데이터 포인트 플롯
for i in range(len(calc_systems)):
    ax.plot(calc_systems[i][1]-calc_systems[i][0], -(calc_systems[i][4]-no3rr_onset),
            mec=calc_systems[i][6], mfc=calc_systems[i][13], mew=0.8, zorder=4,
            marker=calc_systems[i][12], linestyle='',
            label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]-no3rr_onset)))

# 범례 추가
ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=0.5, ncol=3, 
         fancybox=True, shadow=False, fontsize='x-small', handlelength=2)

fig.savefig('NO3RR_1D.pdf', bbox_inches='tight')
plt.close() 