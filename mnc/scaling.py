for dir in import os
import re
import glob
import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from ase.io import read, write


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

# Define the rows and spins
rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}
adsorbates = ['clean', 'oh', 'o', 'ooh']

def main():
    for row_key, metals in rows.items():
        for m, metal in enumerate(metals):
            for ads in adsorbates:
                if ads == 'clean':
                    path = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/{m+2}_{metal}/most_stable/relaxed'
                else:
                    path = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/{m+2}_{metal}/most_stable/{ads}'
                path_DONE = os.path.join(path, 'DONE')
                path_json = os.path.join(path, 'final_with_calculator.json')
                if os.path.exists(path_DONE) and os.path.exists(path_json):
                    atoms = read(path_json)
                    energy = atoms.get_total_energy()
                    df.at[metal, ads] = energy
                else:
                    df.at[metal, ads] = np.nan
    print(df)
                    
if __name__ == '__main__':
    main()
