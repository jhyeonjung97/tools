import os
import re
from statistics import mean
import glob
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
from ase.io import read, write
from matplotlib import rc

# Define the rows and spins
rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
dirs = ['relaxed', 'oh', 'o', 'ooh']
adsorbates = ['clean', 'OH', 'O', 'OOH']
df = pd.DataFrame()

def main():
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals):
            for d, dir in enumerate(dirs):
                ads = adsorbates[d]
                path = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/{m+2}_{metal}/most_stable/{dir}'
                path_DONE = os.path.join(path, 'DONE')
                path_json = os.path.join(path, 'final_with_calculator.json')
                if os.path.exists(path_DONE) and os.path.exists(path_json):
                    atoms = read(path_json)
                    energy = atoms.get_total_energy()
                    df.at[metal, {ads}] = energy
                else:
                    df.at[metal, {ads}] = np.nan

    df['G_'] = df['clean']
    df['G_OH'] = df['OH'] + OH_corr
    df['G_O'] = df['O'] + O_corr
    df['G_OOH'] = df['OOH'] + OOH_corr

                    
if __name__ == '__main__':
    
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
    
    OOH_Cv = 0.077
    OOH_TS = 0.134
    OOH_ZPE = 0.471
    OOH_corr = OOH_Cv - OOH_TS + OOH_ZPE
    
    main()
    
