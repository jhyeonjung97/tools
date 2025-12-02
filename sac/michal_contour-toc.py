import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from matplotlib.colors import ListedColormap

# Global settings for plot
def set_plot_style():
    plt.rcParams.update({
        'ps.usedistiller': 'xpdf',
        'font.size': 9,
        'axes.labelsize': 9,
        'legend.fontsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth': 1.0,
        'axes.linewidth': 0.8,
        'text.usetex': False  # Disable LaTeX rendering
    })

# Function to calculate OOH scaling
def ooh_oh_scaling(doh):
    return 0.84 * doh + 3.14

# Function to calculate overpotential
def overpotential_oer(x, doh):
    dooh = ooh_oh_scaling(doh)
    dg14 = [doh, x, dooh - (x + doh), -dooh + 4.92]
    return max(dg14) - 1.23

# Function for plotting contour and colorbars
def plot_contour(X, Y, Z, x1, x2, y1, y2, fig, ax):
    levels = np.arange(0.2, 1.6, 0.1)
    cmap = ListedColormap([
        '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
        '#ffffe5', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
    ])
    CS = ax.contourf(X, Y, Z, levels, cmap=cmap, extend='max', origin='lower')
    
    # Colorbar
    cbar = plt.colorbar(CS, ticks=list(np.arange(0.2, 1.6, step=0.1)))
    cbar.ax.set_ylabel(r'$\eta_{\sf OER}$ (V)')
    cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

    # Axis settings
    ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize=10)
    ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)', fontsize=10)
    ax.axis([x1, x2, y1, y2])

# Plot data points for systems
def plot_systems(ax, calc_systems):
    for system in calc_systems:
        ax.plot(system[1] - system[0], system[0], 'o', mec=system[5], mfc=system[5], 
                mew=0.8, zorder=4, label=f'{system[4]}: {system[3]:.2f} V')

# Function to save overpotential results into a .tsv file
def save_to_tsv(df, filename='contour_OER_results.tsv'):
    df.to_csv(filename, sep=',', index=False)
    print(f"Data saved to {filename}")

# Main plotting function
def plot_oer_contour():
    set_plot_style()
    
    fig_width_pt = 1.8 * 246.0
    inches_per_pt = 1.0 / 72.27
    golden_mean = (np.sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean
    fig_size = [fig_width, fig_height]
    
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    
    zoom = 0.5
    xcenter, ycenter = 1.45, 0.73
    d1, d2 = 3 * zoom, 4 * zoom
    x1, x2 = xcenter - d1, xcenter + d1
    y1, y2 = ycenter - d2, ycenter + d2
    
    delta = 0.01
    x = np.arange(x1, x2 + delta, delta)
    y = np.arange(y1, y2 + delta, delta)    
    X, Y = np.meshgrid(x, y)
    Z = np.array([[overpotential_oer(i, j) for i in x] for j in y])

    plot_contour(X, Y, Z, x1, x2, y1, y2, fig, ax)
    
    # Load data from 'scaling_relationship.tsv'
    df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figure/scaling_relationship.tsv', sep=',', index_col=0)

    # Add calculated dOOH and overpotential
    df['dG_OOH'] = df['dG_OH'].apply(ooh_oh_scaling)
    df['overpotential'] = df.apply(lambda row: overpotential_oer(row['dG_O'], row['dG_OH']), axis=1)

    # Extract data for plotting (replace dummy data)
    calc_systems = []
    for idx, row in df.iterrows():
        calc_systems.append([row['dG_OH'], row['dG_O'], row['dG_OOH'], row['overpotential'], 
                             f'System {idx}', '#737373', 0.0, 0.08, 1.5, 'o', 'blue', 'red'])
    
    # Plot the systems
    plot_systems(ax, calc_systems)
    
    # Save figure
    fig.savefig('contour_OER.png', bbox_inches='tight')
    fig.clf()
    
    # Save results to a .tsv file
    save_to_tsv(df)

# Run the plot function
plot_oer_contour()
