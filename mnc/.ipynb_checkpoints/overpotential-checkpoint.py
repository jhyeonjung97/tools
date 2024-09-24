import os
import glob
import numpy as np
import pandas as pd
from ase.io import read
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import make_interp_spline

steps = ['H2O->*OH','*OH->*O', '*O->*OOH', '*OOH->O2']
rows = ['3d', '3d', '3d', '3d', '4d', '5d']
groups = ['5', '6', '7', '8', '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
adsorbates = ['clean', 'O', 'OH']

root = '/pscratch/sd/j/jiuy97/6_MNC/figure'
for i in range(6):
    row = rows[i]
    group = groups[i]
    metal = metals[i]
    energies = {}
    for adsorbate in adsorbates:
        tsv_path = os.path.join(root, f'{row}_{group}{metal}_{adsorbate}.tsv')
        energies[adsorbate] = pd.read_csv(tsv_path, sep='\t', index_col=0)
        energies[adsorbate] = energies[adsorbate].head(7)

        if 'LS' in energies[adsorbate].columns:
            energies[adsorbate] = energies[adsorbate].drop(columns=['LS'])
        if 'IS' in energies[adsorbate].columns:
            energies[adsorbate] = energies[adsorbate].drop(columns=['IS'])
        if 'HS' in energies[adsorbate].columns:
            energies[adsorbate] = energies[adsorbate].drop(columns=['HS'])

    
    print(f"Energy data for {metal}:")
    print(energies)