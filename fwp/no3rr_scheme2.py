#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

fig_size = (6, 4)
x1, x2 = -1.0, 1.5
y1, y2 = -3.5, 0.2
delta = 0.01

plt.figure(figsize=fig_size, dpi=300)

xx = np.arange(x1, x2, delta)

plt.plot(xx, -np.maximum(-2*xx, 1.4*xx+0.6), '--', color='black', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx, -np.maximum(-2*xx, 1.6*xx+0.7), '--', color='black', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx, -np.maximum(-2*xx, 1.8*xx+0.8), '--', color='black', lw=0.67, dashes=(3, 1), zorder=2)

plt.plot(xx, -np.maximum(-xx, xx), '--', color='silver', lw=0.67, dashes=(3, 1), zorder=2)
plt.plot(xx, np.ones_like(xx)*0.056, '-', color='black', lw=0.67, zorder=2)

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
    plt.scatter(x, -max(-x, x), color='black', marker=systems[i][3], linewidth=0.67, zorder=3, s=12)
    plt.scatter(x, -max(-2*x, 1.4*x+0.6), facecolor=systems[i][2], edgecolor='black', marker=systems[i][3], linewidth=0.67, zorder=3)
    if -2*x < 1.6*x+0.7:
        plt.scatter(x, -max(-2*x, 1.6*x+0.7), facecolor='white', edgecolor='black', marker=systems[i][3], linewidth=0.67, zorder=3, s=12)
        plt.scatter(x, -max(-2*x, 1.8*x+0.8), facecolor='white', edgecolor='black', marker=systems[i][3], linewidth=0.67, zorder=3, s=12)
    if systems[i][0] == 'CuN4':
        plt.annotate(systems[i][0], xy=(x, -max(-2*x, 1.4*x+0.6)), xytext=(x+0.01, -max(-2*x, 1.4*x+0.6)+0.06), fontsize=8, ha='center')
    elif systems[i][0] in ['RhN4', 'TiN4']:
        plt.annotate(systems[i][0], xy=(x, -max(-2*x, 1.4*x+0.6)), xytext=(x-0.05, -max(-2*x, 1.4*x+0.6)), fontsize=8, ha='right')
    elif systems[i][0] == 'Cu(100)':
        plt.annotate(systems[i][0], xy=(x, -max(-2*x, 1.4*x+0.6)), xytext=(x+0.05, -max(-2*x, 1.4*x+0.6)), fontsize=8)
    else:
        plt.annotate(systems[i][0], xy=(x, -max(-2*x, 1.4*x+0.6)), xytext=(x, -max(-2*x, 1.4*x+0.6)+0.05), fontsize=8)

plt.text(0.5, -0.1, r'E°$_{N_{2}/NH_{3}}$ = 0.05 V vs RHE', color='black', fontsize=8)
# plt.xticks([])
# plt.yticks([])
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlim(x1, x2)
plt.ylim(y1, y2)
plt.xlabel(r'$\Delta G_{\sf H}$ (eV)')
plt.ylabel(r'$U_{\sf NO3RR}$ (V)')
plt.tight_layout()
plt.savefig('no3rr_scheme.png', dpi=300, bbox_inches='tight', transparent=True)
plt.close()