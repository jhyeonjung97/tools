import matplotlib.pyplot as plt
from ase.io import read
from statistics import mean
import glob
import os
import pandas as pd
from scipy.interpolate import make_interp_spline
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

custom_legend = [
    Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#ff7f0e', 
           label='LS (constrained Δz, enforced spin)', markersize=8, linestyle='None'),
    Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#279ff2', 
           label='IS (constrained Δz, enforced spin)', markersize=8, linestyle='None'),
    Line2D([0], [1], marker='s', markerfacecolor='white', markeredgecolor='#9467bd', 
           label='HS (constrained Δz, enforced spin)', markersize=8, linestyle='None'),

    Line2D([0], [0], marker='s', color='#ff7f0e', 
           label='LS (relaxed Δz, free spin)', markersize=8, linestyle='None'),
    Line2D([0], [0], marker='s', color='#279ff2', 
           label='IS (relaxed Δz, free spin)', markersize=8, linestyle='None'),
    Line2D([0], [0], marker='s', color='#9467bd', 
           label='HS (relaxed Δz, free spin)', markersize=8, linestyle='None'),

    Line2D([0], [0], marker='x', color='#ff7f0e', 
           label='LS (constrained Δz, free spin)', markersize=8, linestyle='None'),
    Line2D([0], [0], marker='x', color='#279ff2', 
           label='IS (constrained Δz, free spin)', markersize=8, linestyle='None'),
    Line2D([0], [0], marker='x', color='#9467bd', 
           label='HS (constrained Δz, free spin)', markersize=8, linestyle='None'),
]
fig, ax = plt.subplots()
ax.legend(handles=custom_legend, ncol=3)
ax.axis('off')
png_filename = "custom_legend_markers_only.png"  # Update with your file path
plt.savefig(png_filename, bbox_inches="tight")
print(f"Figure saved as {png_filename}")
plt.close()

fig, ax = plt.subplots()
ax.legend(handles=custom_legend, ncol=1)
ax.axis('off')
png_filename = "custom_legend_markers_only2.png"  # Update with your file path
plt.savefig(png_filename, bbox_inches="tight")
print(f"Figure saved as {png_filename}")
plt.close()