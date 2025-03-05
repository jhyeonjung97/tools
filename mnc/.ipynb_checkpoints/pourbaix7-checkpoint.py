import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

dirs = ["/pscratch/sd/j/jiuy97/6_MNC/pourbaix/1_Fe/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/2_Co/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/3_Mo/"]
main_dirs = ["clean", "mh", "nh", "oh", "o",
             "ohoh", "oh-oh", "ohooh", "oohoh", "oh-ooh", "ooh-oh",
             "ooh", "oho", "oh-o", "o-oh", "oo", "o-o",
             "oooh", "ooho", "o-ooh", "ooh-o", "oohooh", "ooh-ooh",
             "oo", "oo-oo", "oo-ohh", "oo-o", "oo-oh", "oo-ooh",
             "ohh", "ohh-oo", "ohh-ohh", "ohh-o", "ohh-oh", "ohh-ooh",]
sub_dirs = ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"]

e0_pattern = re.compile(r"E0=\s*(-?\.\d+E[+-]?\d+)")
energy_pattern = re.compile(r"energy\s*=\s*(-?\d+\.\d+)")
    
fig_width_pt = 1.8 * 246.0  # LaTeX column width
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # Width in inches
fig_height = fig_width * golden_mean  # Height in inches
fig_size = [fig_width, fig_height]

font_size = 10
tick_font_size = 10

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = 2 * font_size
plt.rcParams['legend.fontsize'] = font_size
plt.rcParams['xtick.labelsize'] = tick_font_size
plt.rcParams['ytick.labelsize'] = tick_font_size
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['lines.linewidth'] = 1.0

# Constants and variables for plotting
Umin, Umax = -0.5, 2.5
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485

pH = np.arange(0, 14, 0.10)
U = np.arange(Umin, Umax, 0.01)
Umax2 = Umax + 0.06 * 14
U2 = np.arange(Umin, Umax2, 0.01)

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

gh = gh2 / 2
go = gh2o - gh2
go2 = 2*gh2o - 2*gh2 + 4.92
print(go, go2/2)
goh = gh2o - gh2 / 2
gooh = 2 * gh2o - 1.5 * gh2

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
dgh2o = dgoh + dgh

color = ['darkgray', ##
         'cornflowerblue', ## 
         'yellowgreen', 
         'teal', 
         'tan', ##
         'salmon', ##
         'forestgreen', 
         'lightsteelblue', ##
         'orange', 
         'gold', 
         'pink', ##
         'plum', ##
         'navy',
         'darkgray', ##
         'cornflowerblue', ## 
         'yellowgreen', 
         'teal', 
         'tan', ##
         'salmon', ##
         'forestgreen', 
         'lightsteelblue', ##
         'orange', 
         'gold', 
         'pink', ##
         'plum', ##
         'navy']
pH2 = np.arange(0, 14.01, 0.01)

figure_path = '/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix'
metal_path = '~/bin/tools/mnc/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
metals = ['Fe', 'Co', 'Mo']

# Function to findthe lowest E0 value in each subdirectory
def find_min_e0(main_dir, sub_dirs):
    min_e0 = None
    for sub_dir in sub_dirs:
        oszicar_path = os.path.join(main_dir, sub_dir, "OSZICAR")
        if os.path.isfile(oszicar_path):
            with open(oszicar_path, 'r') as file:
                for line in file:
                    match = e0_pattern.search(line)
                    if match:
                        e0_value = float(match.group(1))
                        if min_e0 is None or e0_value < min_e0:
                            min_e0 = e0_value
    return min_e0
    
# Function to extract energy from DONE in most_stable or find min_e0 as fallback
def get_energy(main_dir, sub_dirs, metal):
    min_e0 = None
    if metal == 'Fe':
        i=3; j=6; k=2
    elif metal == 'Co':
        i=3; j=7; k=3
    elif metal == 'Mo':
        i=4; j=4; k=5
    if main_dir == 'clean':
        most_stable_dir = f"/pscratch/sd/j/jiuy97/6_MNC/0_clean/{i}d/{j}_{metal}/most_stable/relaxed"
    elif main_dir == 'o':
        most_stable_dir = f"/pscratch/sd/j/jiuy97/6_MNC/1_O/{k}_{metal}/most_stable/relaxed"
    elif main_dir == 'oh':
        most_stable_dir = f"/pscratch/sd/j/jiuy97/6_MNC/2_OH/{k}_{metal}/most_stable/relaxed"
    elif main_dir == 'ooh':
        most_stable_dir = f"/pscratch/sd/j/jiuy97/6_MNC/3_OOH/{k}_{metal}/most_stable/relaxed"
    else:
        most_stable_dir = os.path.join(main_dir, "most_stable")
    done_path = os.path.join(most_stable_dir, "DONE")
    json_path = os.path.join(most_stable_dir, "final_with_calculator.json")
    if os.path.isfile(done_path) and os.path.isfile(json_path):
        atoms = read(json_path)
        min_e0 = atoms.get_potential_energy()
        return min_e0
    # Fallback to finding minimum E0 value
    return find_min_e0(main_dir, sub_dirs)
    if min_e0 is None:
        print(f"Warning: No valid energy found in {main_dir}.")
    return min_e0
    
def addH(x, y):
    return -gh + dgh + 1 * (y + x * const)
    
def addO(x, y):
    return -go + dgo - 2 * (y + x * const)

def addOH(x, y):
    return -goh + dgoh - (y + x * const)

def addOOH(x, y):
    return -gooh + dgooh - 3 * (y + x * const)

def addH2O(x, y):
    return -gh2o + dgh2o

def addO2(x, y):
    return -go2 + dgo2

def dg(i, x, y):
    if surfs[i][0] is None:
        return None
    dg = (surfs[i][0]
          - surfs[1][0]
          + surfs[i][1] * addH(x, y)
          + surfs[i][2] * addO(x, y) 
          + surfs[i][3] * addOH(x, y) 
          + surfs[i][4] * addOOH(x, y)
          + surfs[i][5] * addH2O(x, y)
          + surfs[i][6] * addO2(x, y))
    if i == 0 and surfs[i][1] == 2:
        return dg + bulk_metal
    return dg
        
def overpotential(int1, int2, int3, int4, df, OER, ORR):
    ints = [int1, int2, int3, int4]
    for i, int in enumerate(ints):
        if isinstance(int, tuple):
            if np.isnan(df.loc[int[1], 'E']):
                ints[i] = int[0]
            elif np.isnan(df.loc[int[0], 'E']):
                ints[i] = int[1]
            elif df.loc[int[0], 'E'] < df.loc[int[1], 'E']:
                ints[i] = int[0]
            else:
                ints[i] = int[1]
            # print(int, ints[i])
    int1, int2, int3, int4 = ints
    
    dG12 = df.loc[int2, 'dG'] - df.loc[int1, 'dG']
    dG23 = df.loc[int3, 'dG'] - df.loc[int2, 'dG']
    dG34 = df.loc[int4, 'dG'] - df.loc[int3, 'dG']
    dG41 = 4.92 - dG12 - dG23 - dG34
    if any(np.isnan(value) for value in [dG12, dG23, dG34, dG41]):
        onsetP_oer = np.nan
        overP_oer = np.nan
        onsetP_orr = np.nan
        overP_orr = np.nan
    else:
        onsetP_oer = max(dG12, dG23, dG34, dG41)
        overP_oer = onsetP_oer - 1.23
        onsetP_orr = min(dG12, dG23, dG34, dG41)
        overP_orr = 1.23 - onsetP_orr
        
    OER['int1'].append(int1); OER['int2'].append(int2); OER['int3'].append(int3); OER['int4'].append(int4)
    OER['dg12'].append(dG12); OER['dg23'].append(dG23); OER['dg34'].append(dG34); OER['dg41'].append(dG41)
    OER['overP'].append(overP_oer); OER['onsetP'].append(onsetP_oer)
    
    ORR['int1'].append(int4); ORR['int2'].append(int3); ORR['int3'].append(int2); ORR['int4'].append(int1)
    ORR['dg12'].append(1.23-dG41+1.23); ORR['dg23'].append(1.23-dG34+1.23); ORR['dg34'].append(1.23-dG23+1.23); ORR['dg41'].append(1.23-dG12+1.23)
    ORR['overP'].append(overP_orr); ORR['onsetP'].append(onsetP_orr)
    
for dir in dirs:
    os.chdir(dir)
    basename = os.path.basename(os.path.normpath(dir))
    A, B = basename.split('_', 1)
    bulk_metal = metal_df.at[B, 'energy']

    df = pd.DataFrame()
    OER = {'int1': [], 'int2': [], 'int3': [], 'int4': [], 'dg12': [], 'dg23': [], 'dg34': [], 'dg41': [], 'overP': [], 'onsetP': []}
    ORR = {'int1': [], 'int2': [], 'int3': [], 'int4': [], 'dg12': [], 'dg23': [], 'dg34': [], 'dg41': [], 'overP': [], 'onsetP': []}
    
    json_path = '/pscratch/sd/j/jiuy97/6_MNC/empty/2_/final_with_calculator.json'
    atoms = read(json_path)
    df.loc['vac', 'E'] = atoms.get_potential_energy()
    
    # Iterate through each main directory to extract E0 values and plot
    for main_dir in main_dirs:
        min_e0 = get_energy(main_dir, sub_dirs, B)
    
        if min_e0 is None:
            # print(f"Missing data in directory '{main_dir}' for plotting.")
            df.loc[main_dir, 'E'] = np.nan
            continue
        else:
            df.loc[main_dir, 'E'] = min_e0
    
    df.loc['vac', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [2, 0, 0, 0, 0, 0]
    df.loc['clean', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 0, 0] # [energy, #Hs, #Os, #OHs, #OOHs]
    df.loc['mh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [1, 0, 0, 0, 0, 0]
    df.loc['nh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [1, 0, 0, 0, 0, 0]
    df.loc['o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 0, 0, 0]
    df.loc['oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 0, 0, 0]
    df.loc['ohoh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 2, 0, 0, 0]
    df.loc['oh-oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 2, 0, 0, 0]
    df.loc['ohooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 1, 0, 0]
    df.loc['oohoh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 1, 0, 0]
    df.loc['oh-ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 1, 0, 0]
    df.loc['ooh-oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 1, 0, 0]
    df.loc['ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 1, 0, 0]
    df.loc['oho', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 1, 0, 0, 0]
    df.loc['oh-o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 1, 0, 0, 0]
    df.loc['o-oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 1, 0, 0, 0]
    df.loc['oo', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 2, 0, 0, 0, 0]
    df.loc['o-o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 2, 0, 0, 0, 0]
    df.loc['oooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 1, 0, 0]
    df.loc['ooho', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 1, 0, 0]
    df.loc['o-ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 1, 0, 0]
    df.loc['ooh-o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 1, 0, 0]
    df.loc['oohooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 2, 0, 0]
    df.loc['ooh-ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 2, 0, 0]
    
    df.loc['oo', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 0, 1]
    df.loc['oo-oo', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 0, 2]
    df.loc['oo-ohh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 1, 1]
    df.loc['oo-o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 0, 0, 1]
    df.loc['oo-oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 0, 0, 1]
    df.loc['oo-ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 1, 0, 1]
    
    df.loc['ohh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 1, 0]
    df.loc['ohh-oo', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 1, 1]
    df.loc['ohh-ohh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 0, 2, 0]
    df.loc['ohh-o', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 1, 0, 0, 1, 0]
    df.loc['ohh-oh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 1, 0, 1, 0]
    df.loc['ohh-ooh', ['#H', '#O', '#OH', '#OOH', '#H2O', '#O2']] = [0, 0, 0, 1, 1, 0]
    
    df['G'] = (df['E'] 
               + dgh * df['#H'] 
               + dgo * df['#O'] 
               + dgoh * df['#OH'] 
               + dgooh * df['#OOH'] 
               + dgh2o * df['#H2O'] 
               + dgo2 * df['#O2'])
    df['dG'] = (df['G'] 
                - df.loc['clean', 'E'] 
                - gh * df['#H'] 
                - go * df['#O'] 
                - goh * df['#OH'] 
                - gooh * df['#OOH'] 
                - gh2o * df['#H2O'] 
                - go2 * df['#O2'])
    df.loc['vac', 'dG'] += bulk_metal
    
    if A == '1' and B == 'Fe':
        overpotential('clean', 'oh', 'o', 'ooh', df, OER, ORR)
        overpotential('oh', 'oh-oh', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, OER, ORR)
        overpotential('oh', 'ohoh', 'oho', ('ohooh', 'oohoh'), df, OER, ORR)
        overpotential('oo', 'oo-oh', 'oo-o', 'oo-ooh', df, OER, ORR)
        overpotential('ohh', 'ohh-oh', 'ohh-o', 'ohh-ooh', df, OER, ORR)
    elif A == '2' and B == 'Co':
        overpotential('clean', 'oh', 'o', 'ooh', df, OER, ORR)
        overpotential('oh', 'oh-oh', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, OER, ORR)
        overpotential('oh', 'ohoh', 'oho', ('ohooh', 'oohoh'), df, OER, ORR)
        overpotential('oo', 'oo-oh', 'oo-o', 'oo-ooh', df, OER, ORR)
        overpotential('ohh', 'ohh-oh', 'ohh-o', 'ohh-ooh', df, OER, ORR)
    elif A == '3' and B == 'Mo':
        overpotential('o', 'oho', 'oo', ('oooh', 'ooho'), df, OER, ORR)

    # Define surfaces with extracted E0 values
    surfs = [
        df.loc['vac', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),
        df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #1
        df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #2
        df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #3
        df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #4
        df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #5
        df.loc['o-oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #6
        df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #7
        df.loc['o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #8
        df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #9
        # df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),
        # df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),
        # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),
        # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),
        df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(), #10
        # df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['oo-oo', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['oo-ohh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['oo-o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['oo-oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['oo-ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['ohh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['ohh-oo', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['ohh-ohh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['ohh-o', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        df.loc['ohh-oh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
        # df.loc['ohh-ooh', ['E', '#H', '#O', '#OH', '#OOH', '#H2O', '#O2']].tolist(),  
    ]
        
    surfs = [surf for surf in surfs if not any(pd.isna(x) for x in surf)]
    
    nsurfs = len(surfs)
    lowest_surfaces = []
    
    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs) if dg(k, 0, j) is not None]
        # if -0.01 < j and j < 0.01:
        #     print(j, values)
        lowest_surfaces.append(np.argmin(values))
        
    crossover = []
    uniquesurf = [lowest_surfaces[0]]
    old_value = lowest_surfaces[0]
    crossover.append(Umin)
    
    for j in range(len(U2)):
        if lowest_surfaces[j] != old_value:
            uniquesurf.append(lowest_surfaces[j])
            crossover.append(U2[j])
            old_value = lowest_surfaces[j]
    crossover.append(Umax2)
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([0, 14, Umin, Umax])
    ax.set_xlabel(r'pH', fontsize='large')
    ax.set_ylabel(r'U/V', fontsize='large')
    current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    extraticks = [1.23]
    combined_ticks = sorted(set(current_yticks) | set(extraticks))
    plt.yticks(combined_ticks)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i H2O-%i O2-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4], surfs[k][5], surfs[k][6])
        plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                         facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)
    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    plt.savefig(f'{figure_path}/{A}{B}_pourbaix_full_oo.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_full_oo.png")
    plt.close()

    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([0, 14, Umin, Umax])
    ax.set_xlabel(r'pH', fontsize='large')
    ax.set_ylabel(r'U/V', fontsize='large')
    current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    extraticks = [1.23]
    combined_ticks = sorted(set(current_yticks) | set(extraticks))
    plt.yticks(combined_ticks)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                         facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[k], alpha=0.3, linewidth=5)
    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    plt.savefig(f'{figure_path}/{A}{B}_pourbaix_clean_oo.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_clean_oo.png")
    plt.close()
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([-1.0, 2.5, -600, 200])
    ax.set_xlabel(r'RHE (V)', fontsize='large')
    ax.set_ylabel(r'$\Delta$G (kJ/mol)', fontsize='large')
    xx = np.arange(-1.00, 2.55, 0.01)
    for k in range(nsurfs):
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i H2O-%i O2-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4], surfs[k][5], surfs[k][6])
        dg_value = dg(k, 0, xx)
        if dg_value is not None:
            ax.plot(xx, dg_value * kjmol, '-', lw=1, c=color[k], label=label)
        else:
            print(f"Skipping plot for surface {k} due to missing data.")
    plt.xlim(-1.0, 2.5)
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    plt.savefig(f'{figure_path}/{A}{B}_pourbaix_oo.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_oo.png")
    plt.close()

    df['#H'] = df['#H'].astype(int)
    df['#O'] = df['#O'].astype(int)
    df['#OH'] = df['#OH'].astype(int)
    df['#OOH'] = df['#OOH'].astype(int)
    df['#H2O'] = df['#H2O'].astype(int)
    df['#O2'] = df['#O2'].astype(int)
    df.to_csv(f'{figure_path}/{A}{B}_energies_oo.csv', sep=',') #, index=False)
    df.to_csv(f'{figure_path}/{A}{B}_energies_oo.tsv', sep='\t', float_format='%.2f')
    print(f"Data saved as {A}{B}_energies_oo.png")

    OER_df = pd.DataFrame(OER)
    OER_df.to_csv(f'{figure_path}/{A}{B}_oer_oo.csv', sep=',') #, index=False)
    OER_df.to_csv(f'{figure_path}/{A}{B}_oer_oo.tsv', sep='\t', float_format='%.2f')
    print(f"Data saved as {A}{B}_oer_oo.png")

    ORR_df = pd.DataFrame(ORR)
    ORR_df.to_csv(f'{figure_path}/{A}{B}_orr_oo.csv', sep=',') #, index=False)
    ORR_df.to_csv(f'{figure_path}/{A}{B}_orr_oo.tsv', sep='\t', float_format='%.2f')
    print(f"Data saved as {A}{B}_orr_oo.png")
    