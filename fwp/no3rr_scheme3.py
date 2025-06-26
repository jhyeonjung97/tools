#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

fig_size = (6, 4)
x1, x2 = -1.5, 1.5
y1, y2 = -2.5, 0.2
delta = 0.01

plt.figure(figsize=fig_size, dpi=300)

xx = np.arange(x1, x2, delta)
xx_left = np.arange(x1, -11/17, delta)

# 선 그리기
yy1 = -np.maximum(-2*xx+0.4, 1.4*xx+2.6)
yy2 = -np.maximum(-2*xx+0.4, 1.6*xx+3.0)
yy3 = -np.maximum(-2*xx+0.4, 1.8*xx+3.4)

plt.plot(xx, yy1, '--', color='silver', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx, yy2, '--', color='silver', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx, yy3, '--', color='silver', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx_left, -(-2*xx_left+0.4), '--', color='black', lw=0.67, dashes=(3, 1), zorder=3)

# 선들 사이의 공간을 색칠
plt.fill_between(xx, yy1, yy3, alpha=0.2, color='silver', zorder=1)

plt.plot(xx, -np.maximum(-xx, xx), '--', color='silver', lw=0.67, dashes=(3, 1), zorder=2)
# plt.plot(xx, np.ones_like(xx)*0.05, '-', color='black', lw=0.67, zorder=2)

systems = []
hbe_corr = 0.27
systems.append(['Cu(111)', 0.01+hbe_corr, 'gold', 'o'])
systems.append(['Cu(100)', 0.09+hbe_corr, 'gold', '^'])
systems.append(['Ni(111)', -0.68+hbe_corr, 'green', 'o'])
# systems.append(['Ru(0001)', -0.88+hbe_corr, 'mediumorchid', 'v'])
systems.append(['Pt(111)', -0.18+hbe_corr, 'dodgerblue', 'o'])
systems.append(['Rh(111)', -0.32+hbe_corr, 'salmon', 'o'])
systems.append(['Ag(111)', 0.38+hbe_corr, 'gray', 'o'])
systems.append(['TiN4', -0.7, 'cyan', 's'])
# systems.append(['RuN4', -0.8, 'magenta', 's'])
systems.append(['RhN4', -0.85, 'red', 's'])
systems.append(['CuN4', 1.4, 'yellow', 's'])
systems.append(['NiN4', 1.2, 'lime', 's'])

for i in range(len(systems)):
    x = systems[i][1]
    plt.scatter(x, -max(-x, x), color='silver', marker=systems[i][3], linewidth=0.67, zorder=4, s=12)
    plt.scatter(x, -max(-2*x+0.4, 1.4*x+2.6), facecolor=systems[i][2], edgecolor='black', marker=systems[i][3], linewidth=0.67, zorder=4)
    if -x+0.2 < 0.9*x+1.5:
        plt.scatter(x, -max(-2*x+0.4, 1.6*x+3.0), facecolor='white', edgecolor='silver', marker=systems[i][3], linewidth=0.67, zorder=4, s=12)
        plt.scatter(x, -max(-2*x+0.4, 1.8*x+3.4), facecolor='white', edgecolor='silver', marker=systems[i][3], linewidth=0.67, zorder=4, s=12)
    # if systems[i][0] == 'CuN4':
    #     plt.annotate(systems[i][0], xy=(x, -max(-x+0.2, 0.7*x+1.3)), xytext=(x+0.01, -max(-x+0.2, 0.7*x+1.3)+0.06), fontsize=8, ha='center')
    # elif systems[i][0] in ['RhN4', 'TiN4']:
    #     plt.annotate(systems[i][0], xy=(x, -max(-x+0.2, 0.7*x+1.3)), xytext=(x-0.05, -max(-x+0.2, 0.7*x+1.3)), fontsize=8, ha='right')
    # elif systems[i][0] == 'Cu(100)':
    #     plt.annotate(systems[i][0], xy=(x, -max(-x+0.2, 0.7*x+1.3)), xytext=(x+0.05, -max(-x+0.2, 0.7*x+1.3)), fontsize=8)
    # else:
    #     plt.annotate(systems[i][0], xy=(x, -max(-x+0.2, 0.7*x+1.3)), xytext=(x, -max(-x+0.2, 0.7*x+1.3)+0.05), fontsize=8)

plt.text(-1.4, 0.1, r'E°$_{NO_{3}^{-}/NH_{3}}$ = 0.70 V vs RHE', color='black', fontsize=8)
# plt.xticks([])
# plt.yticks([])
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlim(x1, x2)
plt.ylim(y1*2, y2*2)
plt.xlabel(r'$\Delta G_{\sf H}$ (eV)')
plt.ylabel(r'U$_{\text{HER}}$, U$_{\text{NO}_3\text{RR}}$ (V)')

# marker_fcc1 = Line2D([0], [0], color='black', marker='o', markeredgewidth=0.67, markersize=3, linestyle='None')
# marker_fcc2 = Line2D([0], [0], markerfacecolor='white', markeredgecolor='black', marker='o', markeredgewidth=0.67, markersize=3, linestyle='None')
# marker_mnc1 = Line2D([0], [0], color='black', marker='s', markeredgewidth=0.67, markersize=3, linestyle='None')
# marker_mnc2 = Line2D([0], [0], markerfacecolor='white', markeredgecolor='black', marker='s', markeredgewidth=0.67, markersize=3, linestyle='None')
# plt.legend([(marker_fcc1, marker_fcc2), (marker_mnc1, marker_mnc2)], 
#             ['metal surfaces', 'M-N-C single atom'], 
#             handler_map={tuple: HandlerTuple(ndivide=None)}, 
#             loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig('no3rr_scheme.png', dpi=300, bbox_inches='tight') #, transparent=True)
plt.close()