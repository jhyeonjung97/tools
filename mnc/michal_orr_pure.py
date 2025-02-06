import csv
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.markers import MarkerStyle

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
x1, x2 = xcenter - 1.4, xcenter + 1.8 # 3.2
y1, y2 = ycenter - 1.6, ycenter + 1.6 # 3.2

ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize=10)
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)', fontsize=10)

# Define functions for overpotential calculations
def ooh_oh_scaling(doh):
    return doh + 3.2

def orr_step(i):
    steps = ['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']
    return steps[i]

def overpotential_orr(doh, do, dooh):
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def overpotential_orr_full(doh, do, dooh):
    dg14 = [-4.92 + dooh, -dooh + do, -do + doh, -doh]
    m = max(dg14)
    return [round(m + 1.23, 2), round(-m, 2), orr_step(dg14.index(m))]

def overpotential_orr_for_contour(doh, dooh):
    do = 1.469 * doh + 1.253
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

# Read data from the TSV file
df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figure/scaling_relationship.tsv', sep='\t', index_col=0)

# Extract values from the dataframe
doh_values = df['dG_OH']
do_values = df['dG_O']

# Add `dG_OOH` and calculate overpotential for each entry
df['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
df['overpotential'] = df.apply(lambda row: overpotential_orr(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

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

# Plot data points from the TSV file with their calculated overpotentials
for idx, row in df.iterrows():
    ax.plot(row['dG_OH'], row['dG_OOH'], 'o', label=f'{row.name}: {row["overpotential"]:.2f} V')

# Add scaling line
ax.plot(x, x+3.2, '--', lw=1, dashes=(3, 1), c='black')
ax.text(1.1, 2.4, r'$\Delta$G$_{\sf OOH}$=$\Delta$G$_{\sf OH}$+3.2 eV', color='black', fontsize=10)
ax.legend(bbox_to_anchor=(0.5, 1.1), loc='center', borderaxespad=0.0, ncol=3, fancybox=True, shadow=False, fontsize='x-small', handlelength=2)
fig.savefig('contour_ORR.png', bbox_inches='tight')
print("Figure saved as contour_ORR.png")
fig.clf()

# CSV writing for overpotential results
with open('contour_ORR.csv', 'w', newline='') as myfile:
    fieldnames = ['Surface name', 'dOH', 'dO', 'dOOH', 'overpotential', 'onset potential', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_orr_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'Surface name': row.name, 
            'dOH': row['dG_OH'], 'dO': row['dG_O'], 'dOOH': row['dG_OOH'], 
            'overpotential': recalculated_over[0],
            'onset potential': recalculated_over[1], 
            'PLS': recalculated_over[2]
        })

# TSV writing for overpotential results
with open('contour_ORR.tsv', 'w', newline='') as myfile:
    fieldnames = ['Surf.', 'dOH', 'dO', 'dOOH', 'overP', 'onsetP', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames, delimiter='\t')  # Change delimiter to '\t'
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_orr_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'dOH': round(row['dG_OH'], 2),
            'dO': round(row['dG_O'], 2),
            'dOOH': round(row['dG_OOH'], 2),
            'overP': round(recalculated_over[0], 2),
            'onsetP': round(recalculated_over[1], 2),
            'PLS': recalculated_over[2]
        })
