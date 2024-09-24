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
spins = ['LS', 'IS', 'HS']
colors = ['#ff7f0e', '#279ff2', '#9467bd']

OH_Cv = 0.042
OH_TS = 0.066
OH_ZPE = 0.376
OH_corr = OH_Cv - OH_TS + OH_ZPE

O_Cv = 0.034
O_TS = 0.060
O_ZPE = 0.064
O_corr = O_Cv - O_TS + O_ZPE

OOH_Cv = 0.077
OOH_TS = 0.134
OOH_ZPE = 0.471
OOH_corr = OOH_Cv - OOH_TS + OOH_ZPE

root = '/pscratch/sd/j/jiuy97/6_MNC/figure'
for m, metals in enumerate(metals):
    row = rows[m]
    group = groups[m]
    energies = {}
    gibbs_energies = {}
    for adsorbate in adsorbates:
        tsv_path = os.path.join(root, f'{row}_{group}{metal}_{adsorbate}.tsv')
        energies[adsorbate] = pd.read_csv(tsv_path, sep='\t', index_col=0)
        energies[adsorbate] = energies[adsorbate].head(7)
        
        for spin in spins:
            if spin in energies[adsorbate].columns:
                energies[adsorbate] = energies[adsorbate].drop(columns=[spin])
                
    for ms_spin0 in energies['clean'].columns:
        for ms_spin1 in energies['OH'].columns:
            gibbs_energies[f'{ms_spin0}->{ms_spin1}'] = energies['clean'][ms_spin0] - energies['OH'][ms_spin1]
    
    if metal == 'Mn':
        print(energies)