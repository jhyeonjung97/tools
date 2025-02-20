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

a, b = 0.87, 3.08
c, d = 1.79, 0.62

rows = ['3d', '3d', '3d', '3d'] #, '4d', '5d']
groups = ['5', '6', '7', '8'] #, '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni'] #, 'Mo', 'W']

# Figure and font settings
fig = plt.figure(figsize=(4,4), dpi=300)
font_size = 9
plt.rcParams.update({
    'ps.usedistiller': 'xpdf',
    'font.size': font_size,
    'axes.labelsize': font_size,
    'legend.fontsize': font_size,
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

# Define functions for overpotential calculations
def ooh_oh_scaling(doh):
    return a * doh + b

def orr_step(i):
    steps = ['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']
    return steps[i]

def overpotential_orr(doh, do, dooh):
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def overpotential_orr_for_contour(doh, dooh):
    do = c * doh + d
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def overpotential_orr_full(doh, do, dooh):
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    m = max(dg14)
    return [round(m + 1.23, 2), round(-m, 2), orr_step(dg14.index(m))]
    
# Read data from the TSV file
save_path='/pscratch/sd/j/jiuy97/6_MNC/figures/contour'
df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figures/contour/scaling_relationship.tsv', sep='\t', header=0, index_col=0)

df['overpotential'] = df.apply(lambda row: overpotential_orr(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Prepare separate data for each metal
dfs = {}
for m, metal in enumerate(metals):
    row = rows[m]
    group = groups[m]
    dfs[metal] = pd.read_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/contour/{row}_{group}{metal}_gibbs.tsv', sep='\t', header=0, index_col=0)
    dfs[metal]['overpotential'] = dfs[metal].apply(lambda row: overpotential_orr(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

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

# Plot data points from the TSV file with their calculated overpotentials
markers = ['o', 's', 'd', '^', 'v', '*']  # Different markers for metals
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
    ax.scatter(row.dG_OH, row.dG_OOH, 
               # label=f'{row.Index}: {row.overpotential:.2f} V',               
               s=24, marker='X', 
               linewidths=0.5,
               facecolors=colors[row_num-1],
               edgecolors='black',
               zorder=10)

# Plot the metal-specific data points with colormaps
for m, metal in enumerate(metals):
    for row_num, row in enumerate(dfs[metal].itertuples(), 1):
        ax.scatter(row.dG_OH, row.dG_OOH, 
                   s=24, marker='o', 
                   linewidths=0.5,
                   facecolors=color_ranges[m][row_num-1],
                   edgecolors='black',
                   zorder=9)
for row_num, row in enumerate(df.itertuples(), 1):
    if row_num < 5:
        ax.scatter([], [], 
                   label=f'{row.Index}: {row.overpotential:.2f} V',               
                   s=24, marker='X', 
                   linewidths=0.5,
                   facecolor=colors[row_num-1],
                   edgecolor='black')
ax.scatter([], [], label='relaxed Δz', s=24, marker='X', linewidths=0.5, facecolor='black', edgecolor='black')
ax.scatter([], [], label='fixed Δz', s=24, marker='o', linewidths=0.5, facecolor='black', edgecolor='black')
ax.set_xticks([])
ax.set_yticks([])
ax.plot(x, a * x + b, '--', lw=1, dashes=(3, 1), c='black')
fig.savefig(os.path.join(save_path, 'contour_toc.png'), bbox_inches='tight')
print("Figure saved as contour_toc.png")
fig.clf()