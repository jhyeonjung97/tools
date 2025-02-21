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
# xcenter, ycenter = 1.60, 0.81
# x1, x2 = xcenter - 1.2, xcenter + 1.2 # 3.2
# y1, y2 = ycenter - 2.0, ycenter + 2.0 # 3.6
x1, x2 = 0.2, 2.7
y1, y2 = -1.1, 2.7

ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf O}$ - $\Delta$G$_{\sf OH}$(eV)', fontsize='large')
ax.set_ylabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize='large')

# Define functions for overpotential calculations
def ooh_oh_scaling(doh):
    return a * doh + b

def oer_step(i):
    steps = ['H2O->OH*', 'OH*->O*', 'O*->OOH*', 'OOH*->O2']
    return steps[i]

def overpotential_oer(doh, do, dooh):
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    return max(dg14) - 1.23
    
def overpotential_oer_for_contour(do_doh, doh):
    do = do_doh + doh
    dooh = ooh_oh_scaling(doh)
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    return max(dg14) - 1.23
    
def overpotential_oer_full(doh, do, dooh):
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    m = max(dg14)
    return [m - 1.23, -m, oer_step(dg14.index(m))]
    
# Read data from the TSV file
save_path='/pscratch/sd/j/jiuy97/6_MNC/figures/contour'
df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figures/contour/scaling_relationship.csv', sep=',', index_col=0)

# Extract values from the dataframe
# doh_values = df['dG_OH']
# do_values = df['dG_O']
# dooh_values = df['dG_OOH']
# df['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
df['overpotential'] = df.apply(lambda row: overpotential_oer(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Prepare separate data for each metal
dfs = {}
for m, metal in enumerate(metals):
    row = rows[m]
    group = groups[m]
    dfs[metal] = pd.read_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/contour/{row}_{group}{metal}_gibbs.csv', sep=',', header=0, index_col=0)
    # doh_values = dfs[metal]['dG_OH']
    # do_values = dfs[metal]['dG_O']
    # dooh_values = dfs[metal]['dG_OOH']
    # dfs[metal]['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
    dfs[metal]['overpotential'] = dfs[metal].apply(lambda row: overpotential_oer(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Generate data for contour plot
delta = 0.01
x = np.arange(x1, x2 + delta, delta)
y = np.arange(y1, y2 + delta, delta)
X, Y = np.meshgrid(x, y)
Z = np.array([[overpotential_oer_for_contour(i, j) for i in x] for j in y])

# Plot contour
levels = np.arange(0.2, 1.6, 0.1)
CS = plt.contourf(X, Y, Z, levels, cmap=ListedColormap([
    '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
]), extend='max', origin='lower')
cbar = plt.colorbar(CS, ticks=np.arange(0.2, 1.6, 0.1))
cbar.ax.set_ylabel(r'$\eta_{\sf OER}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

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
    ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
               # label=f'{row.Index}: {row.overpotential:.2f} V',               
               s=24, marker='X', 
               linewidths=0.5,
               facecolor=colors[row_num-1],
               edgecolor='black',
               zorder=10)

# Plot the metal-specific data points with colormaps
for m, metal in enumerate(metals):
    for row_num, row in enumerate(dfs[metal].itertuples(), 1):
        ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
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

ax.plot(x, x + 3.2, '--', lw=1, dashes=(3, 1), c='black')
ax.text(1.1, 2.3, r'$\Delta$G$_{\sf OOH}$=$\Delta$G$_{\sf OH}$+3.2 eV', color='black', fontsize=10)

ax.legend(bbox_to_anchor=(0.5, 1.1), loc='center', borderaxespad=0.5,
          ncol=3, columnspacing=1.0, handletextpad=0.4,
          fancybox=True, shadow=False, fontsize='small', handlelength=2)
# ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1, 5.5, 0.2), borderaxespad=0.5,
#           ncol=2, columnspacing=6.0, handletextpad=0.1,
#           fancybox=True, shadow=False, fontsize='small', handlelength=2)
fig.savefig(os.path.join(save_path, 'contour_OER.png'), bbox_inches='tight')
print("Figure saved as contour_OER.png")
# fig.savefig(os.path.join(save_path, 'contour_OER_label.png'), bbox_inches='tight')
# print("Figure saved as contour_OER_label.png")

fig.clf()

# Water and hydrogen properties
water_E = -14.23797429
water_Cv = 0.103
water_TS = 0.675
water_ZPE = 0.558
water_G = water_E + water_Cv - water_TS + water_ZPE

hydrogen2_E = -6.77412273
hydrogen2_Cv = 0.0905
hydrogen2_TS = 0.408
hydrogen2_ZPE = 0.273
hydrogen2_G = hydrogen2_E + hydrogen2_Cv - hydrogen2_TS + hydrogen2_ZPE
hydrogen_G = hydrogen2_G / 2

oxygen_G = water_G - hydrogen2_G
hydroxide_G = water_G - hydrogen_G

# Correction values
OH_Cv = 0.042
OH_TS = 0.066
OH_ZPE = 0.376
OH_corr = OH_Cv - OH_TS + OH_ZPE

O_Cv = 0.034
O_TS = 0.060
O_ZPE = 0.064
O_corr = O_Cv - O_TS + O_ZPE

E_ = -279.25
E_OH = -289.57160323
E_O = -284.14645182

G_ = E_
G_OH = E_OH + OH_corr
G_O = E_O + O_corr

dG_OH = G_OH - G_ - hydroxide_G
dG_O = G_O - G_ - oxygen_G

for m, metal in enumerate(metals):
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([x1, x2, y1, y2])
    ax.set_xlabel(r'$\Delta$G$_{\sf O}$ - $\Delta$G$_{\sf OH}$(eV)', fontsize='large')
    ax.set_ylabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize='large')
    CS = plt.contourf(X, Y, Z, levels, cmap=ListedColormap([
        '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
        '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
    ]), extend='max', origin='lower')
    # cbar = plt.colorbar(CS, ticks=np.arange(0.3, 1.6, 0.1))
    # cbar.ax.set_ylabel(r'$\eta_{\sf OER}$ (V)')
    # cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')
    row = df.loc[metal]
    ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
               # label=f'{row.name}: {row.overpotential:.2f} V',
               s=36, marker='X', 
               linewidths=0.5,
               facecolor=colors[m],
               edgecolor='black',
               zorder=10)
    for row_num, row in enumerate(dfs[metal].itertuples(), 1):
        ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
                   s=36, marker='o',
                   linewidths=0.5,
                   facecolors=color_ranges[m][row_num-1],
                   edgecolors='black',
                   zorder=9)
    # if metal == 'Co':
    #     ax.scatter(dG_O - dG_OH, dG_OH, 
    #                s=36, marker='o',
    #                linewidths=1.0,
    #                facecolors='white',
    #                edgecolors='black',
    #                zorder=9)
    norm = mcolors.Normalize(vmin=-0.1, vmax=1.3)
    sm = cm.ScalarMappable(norm=norm, cmap=ListedColormap(color_ranges[m]))
    sm.set_array([])
    cbar2 = plt.colorbar(sm, ax=ax, ticks=np.arange(0.0, 1.4, 0.2), extend='max')
    cbar2.ax.set_ylabel(r'$\Delta z$ (Å)')
    cbar2.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')
    fig.savefig(os.path.join(save_path, f"contour_OER_{m+1}{metal}.png"), bbox_inches='tight')
    print(f"Figure saved as contour_OER_{m+1}{metal}.png")
    fig.clf()

# CSV writing for overpotential results
with open(os.path.join(save_path, 'contour_OER.csv'), 'w', newline='') as myfile:
    fieldnames = ['Surface name', 'dOH', 'dO', 'dOOH', 'overpotential', 'onset potential', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'Surface name': row.name, 
            'dOH': row['dG_OH'], 'dO': row['dG_O'], 'dOOH': row['dG_OOH'], 
            'overpotential': recalculated_over[0],
            'onset potential': recalculated_over[1], 
            'PLS': recalculated_over[2]
        })

# TSV writing for overpotential results
with open(os.path.join(save_path, 'contour_OER.csv'), 'w', newline='') as myfile:
    fieldnames = ['Surf.', 'dOH', 'dO', 'dO*', 'diff', 'dOOH', 'overP', 'onsetP', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames) #, delimiter='\t')  # Change delimiter to '\t'
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'Surf.': row.name, 
            'dOH': row['dG_OH'],
            'dO': row['dG_O'],
            'dO*': c*row['dG_OH']+d,
            'diff': c*row['dG_OH']+d-row['dG_O'],
            'dOOH': row['dG_OOH'],
            'overP': recalculated_over[0],
            'onsetP': recalculated_over[1],
            'PLS': recalculated_over[2]
        })

# Write results for each metal
for m, metal in enumerate(metals):
    with open(os.path.join(save_path, f'contour_OER_{m+1}{metal}.csv'), 'w', newline='') as myfile:
        fieldnames = ['Surf.', 'dOH', 'dO', 'dO*', 'diff', 'dOOH', 'overP', 'onsetP', 'PLS']
        writer = csv.DictWriter(myfile, fieldnames=fieldnames) #, delimiter='\t')
        writer.writeheader()
        for idx, row in dfs[metal].iterrows():
            recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
            writer.writerow({
                'Surf.': row.name, 
                'dOH': row['dG_OH'],
                'dO': row['dG_O'],
                'dO*': c*row['dG_OH']+d,
                'diff': c*row['dG_OH']+d-row['dG_O'],
                'dOOH': row['dG_OOH'],
                'overP': recalculated_over[0],
                'onsetP': recalculated_over[1],
                'PLS': recalculated_over[2]
            })
        # writer.writerow({
        #     'Surf.': row.name, 
        #     'dOH': round(row['dG_OH'], 2),
        #     'dO': round(row['dG_O'], 2),
        #     'dO*': round(c*row['dG_OH']+d, 2),
        #     'diff': round(c*row['dG_OH']+d-row['dG_O'], 2),
        #     'dOOH': round(row['dG_OOH'], 2),
        #     'overP': round(recalculated_over[0], 2),
        #     'onsetP': round(recalculated_over[1], 2),
        #     'PLS': recalculated_over[2]
        # })