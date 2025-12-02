import os
import csv
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from matplotlib.markers import MarkerStyle

a, b = 0.87, 3.08 # oh ooh
c, d = 1.79, 0.62 # oh o
e, f = 0.51, -0.78 # oh o2

rows = ['3d', '3d', '3d', '3d'] #, '4d', '5d']
groups = ['5', '6', '7', '8'] #, '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni'] #, 'Mo', 'W']

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
go2 = 2*gh2o - 2*gh2 + 4.92

# print(go2)

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo
dgo2 = dgooh - dgh

# Figure and font settings
fig_width_pt = 1.8 * 246.0
inches_per_pt = 1.0 / 72.27
golden_mean = (np.sqrt(5) - 1.0) / 2.0
fig_width = fig_width_pt * inches_per_pt
fig_height = fig_width * golden_mean
fig_size = [fig_width, fig_height]
fig = plt.figure(figsize=fig_size, dpi=300)

font_size = 9
tick_font_size = 8
plt.rcParams.update({
    'ps.usedistiller': 'xpdf',
    'font.size': font_size,
    'axes.labelsize': font_size,
    'legend.fontsize': font_size,
    'xtick.labelsize': tick_font_size,
    'ytick.labelsize': tick_font_size,
    'lines.linewidth': 1.0
})

def setfont(font='cmss', unicode=True):
    """
    Set Matplotlib rcParams to use LaTeX for font rendering.
    """
    font = font.lower().replace(" ", "")
    font = {'family': 'sans-serif', 'serif': ['cmss']}
    preamble = r"""\usepackage{color} \usepackage[tx]{sfmath}"""
    plt.rc('font', **font)
    plt.rcParams['text.latex.preamble'] = preamble

setfont()

# Plot settings
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
xcenter, ycenter = 1.0, 3.7
x1, x2 = xcenter - 1.4, xcenter + 1.2
y1, y2 = ycenter - 1.6, ycenter + 1.6

ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize='large')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)', fontsize='large')

# Define functions for overpotential calculations (using op04 method from orr_overpotential.py)
def orr_step(i):
    steps = ['O2->OO*', 'OO*->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']
    return steps[i]

def overpotential_orr_op04(doh, do, dooh, do2):
    dg04 = [do2, dooh-4.92-do2+1.23, do-dooh+1.23, doh-do+1.23, -doh+1.23]
    return max(dg04) 

def overpotential_orr_for_contour(doh, dooh):
    do = c * doh + d
    do2 = e * doh + f
    return overpotential_orr_op04(doh, do, dooh, do2)
    
# Read data from the orr_overpotential.py generated CSV file
df = pd.read_csv('orr_results.csv', sep=',', header=0, index_col=0)

# Calculate overpotential using op04 method
df['overpotential'] = df.apply(lambda row: overpotential_orr_op04(row['goh'], row['go'], row['gooh'], row['dg0']), axis=1)

# Prepare separate data for each metal (using the same data for now)
dfs = {}
for m, metal in enumerate(metals):
    dfs[metal] = df[df.index == metal].copy()
    if len(dfs[metal]) > 0:
        dfs[metal]['overpotential'] = dfs[metal].apply(lambda row: overpotential_orr_op04(row['goh'], row['go'], row['gooh'], row['dg0']), axis=1)

# Generate data for contour plot
delta = 0.01
x = np.arange(x1, x2 + delta, delta)
y = np.arange(y1, y2 + delta, delta)
X, Y = np.meshgrid(x, y)
Z = np.array([[overpotential_orr_for_contour(i, j) for i in x] for j in y])

# Plot contour
levels = np.arange(0.1, 1.6, 0.1)
CS = plt.contourf(X, Y, Z, levels, cmap=ListedColormap([
    '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
]), extend='max', origin='lower')
cbar = plt.colorbar(CS, ticks=np.arange(0.1, 1.6, 0.1))
cbar.ax.set_ylabel(r'$\eta_{\sf ORR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# Plot data points from the CSV file with their calculated overpotentials
colors = ['blue', 'green', 'orange', 'red', 'purple', 'grey']
color_ranges = [
    plt.cm.Blues(np.linspace(0.3, 0.9, 7)),
    plt.cm.Greens(np.linspace(0.3, 0.9, 7)),
    plt.cm.Oranges(np.linspace(0.3, 0.9, 7)),
    plt.cm.Reds(np.linspace(0.3, 0.9, 7)),
    plt.cm.Purples(np.linspace(0.3, 0.9, 7)),
    plt.cm.Greys(np.linspace(0.3, 0.9, 7)),
    ]

# Plot the general dataset points
for row_num, row in enumerate(df.itertuples(), 1):
    ax.scatter(row.goh, row.gooh, 
               # label=f'{row.Index}: {row.overpotential:.2f} V',               
               s=24, marker='X', 
               linewidths=0.5,
               facecolors=colors[row_num-1],
               edgecolors='black',
               zorder=10)

for row_num, row in enumerate(df.itertuples(), 1):
    if row_num < 5:
        ax.scatter([], [], 
                   label=f'{row.Index}: {row.overpotential:.2f} V',               
                   s=24, marker='X', 
                   linewidths=0.5,
                   facecolor=colors[row_num-1],
                   edgecolor='black')
ax.scatter([], [], label='relaxed Δz', s=24, marker='X', linewidths=0.5, facecolor='black', edgecolor='black')

ax.plot(x, a * x + b, '--', lw=1, dashes=(3, 1), c='black')
ax.text(1.4, 2.7, r'$\Delta$G$_{\sf OOH}$=', color=(0.8,0.8,0.8), fontsize=10)
ax.text(1.4, 2.5, rf'{a}*$\Delta$G$_{{\sf OH}}$', color=(0.8,0.8,0.8), fontsize=10)
ax.text(1.4, 2.3, rf'+{b} eV', color=(0.8,0.8,0.8), fontsize=10)
ax.legend(bbox_to_anchor=(0.5, 1.1), loc='center', borderaxespad=0.5,
          ncol=3, columnspacing=1.0, handletextpad=0.4,
          fancybox=True, shadow=False, fontsize='small', handlelength=2)
fig.savefig(os.path.join('contour_ORR_o2.png'), bbox_inches='tight')
print("Figure saved as contour_ORR_o2.png")