import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# dirs = glob.glob("/pscratch/sd/j/jiuy97/6_MNC/pourbaix/*_*/")
dirs = ["/pscratch/sd/j/jiuy97/6_MNC/pourbaix/1_Fe/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/2_Co/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/3_Mo/"]
main_dirs = ["clean", "mh", "nh", "oh", "o", 
             "ohoh", "oh-oh", "ohooh", "oohoh", "oh-ooh", "ooh-oh",
             "ooh", "oho", "oh-o", "o-oh", "oo", "o-o",
             "oooh", "ooho", "o-ooh", "ooh-o", "oohooh", "ooh-ooh"]
# main_dirs = ["clean", "h", "oh", "o", 
#              "ohoh", "oh-oh", "ohooh", "oohoh", "ooh-oh", # "oh-ooh"
#              "ooh", "oho", "o-oh", "o-o", "oo", # "oh-o"
#              "oooh", "ooho", "ooh-o", "oohooh", "ooh-ooh"] # "o-ooh"
sub_dirs = ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"]

# Regular expression to match E0 values in scientific notation
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

zpeh2o = 0.560
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2

gh = gh2 / 2
go = gh2o - gh2
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
         'navy']
pH2 = np.arange(0, 14.01, 0.01)

metal_path = '/pscratch/sd/j/jiuy97/6_MNC/gas/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)

elements_data = {
    "Ti": {"electrode_potential": -1.63, "cation_charge": 2},  # Ti^2+ + 2e- → Ti (Ti^2+ prioritized over Ti^4+)
    "V":  {"electrode_potential": -1.18, "cation_charge": 2},  # V^2+ + 2e- → V (V^2+ prioritized over V^3+)
    "Cr": {"electrode_potential": -0.91, "cation_charge": 2},  # Cr^2+ + 2e- → Cr (Cr^2+ prioritized over Cr^3+)
    "Mn": {"electrode_potential": -1.18, "cation_charge": 2},  # Mn^2+ + 2e- → Mn
    "Fe": {"electrode_potential": -0.44, "cation_charge": 2},  # Fe^2+ + 2e- → Fe (Fe^2+ prioritized over Fe^3+)
    "Co": {"electrode_potential": -0.28, "cation_charge": 2},  # Co^2+ + 2e- → Co (Co^2+ prioritized over Co^3+)
    "Ni": {"electrode_potential": -0.25, "cation_charge": 2},  # Ni^2+ + 2e- → Ni
    "Cu": {"electrode_potential": 0.34,  "cation_charge": 2},  # Cu^2+ + 2e- → Cu

    "Zr": {"electrode_potential": -1.45, "cation_charge": 2},  # Zr^2+ + 2e- → Zr (Zr^2+ prioritized over Zr^4+)
    "Nb": {"electrode_potential": -1.00, "cation_charge": 3},  # Nb^3+ + 3e- → Nb (Nb^3+ prioritized over Nb^5+)
    "Mo": {"electrode_potential": -0.20, "cation_charge": 3},  # Mo^3+ + 3e- → Mo (Mo^3+ prioritized over Mo^6+)
    "Tc": {"electrode_potential": 0.74,  "cation_charge": 7},  # Tc^7+ + 7e- → Tc (No lower oxidation state available)
    "Ru": {"electrode_potential": 0.45,  "cation_charge": 2},  # Ru^2+ + 2e- → Ru (Ru^2+ prioritized over Ru^3+)
    "Rh": {"electrode_potential": 0.76,  "cation_charge": 2},  # Rh^2+ + 2e- → Rh (Rh^2+ prioritized over Rh^3+)
    "Pd": {"electrode_potential": 0.95,  "cation_charge": 2},  # Pd^2+ + 2e- → Pd

    "Hf": {"electrode_potential": -1.55, "cation_charge": 4},  # Hf^4+ + 4e- → Hf (No lower oxidation state available)
    "Ta": {"electrode_potential": -1.10, "cation_charge": 3},  # Ta^3+ + 3e- → Ta (Ta^3+ prioritized over Ta^5+)
    "W":  {"electrode_potential": -0.11, "cation_charge": 4},  # W^4+ + 4e- → W (W^4+ prioritized over W^6+)
    "Re": {"electrode_potential": 0.30,  "cation_charge": 4},  # Re^4+ + 4e- → Re (Re^4+ prioritized over Re^7+)
    "Os": {"electrode_potential": 0.40,  "cation_charge": 2},  # Os^2+ + 2e- → Os (Os^2+ prioritized over Os^4+)
    "Ir": {"electrode_potential": 1.16,  "cation_charge": 2},  # Ir^2+ + 2e- → Ir (Ir^2+ prioritized over Ir^3+)
    "Pt": {"electrode_potential": 1.20,  "cation_charge": 2},  # Pt^2+ + 2e- → Pt
}

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
def get_energy(main_dir, sub_dirs):
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

def dg(i, x, y):
    if surfs[i][0] is None:
        return None
    elif i == 0 and surfs[i][1] == 2:
        return (surfs[i][0] 
                - surfs[1][0] 
                + cation - charge * (y + x * const)
                + surfs[i][1] * addH(x, y) 
                + surfs[i][2] * addO(x, y) 
                + surfs[i][3] * addOH(x, y) 
                + surfs[i][4] * addOOH(x, y))
    return (surfs[i][0] 
            - surfs[1][0] 
            + surfs[i][1] * addH(x, y) 
            + surfs[i][2] * addO(x, y) 
            + surfs[i][3] * addOH(x, y) 
            + surfs[i][4] * addOOH(x, y))
    
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
    ORR['dg12'].append(1.23-dG34); ORR['dg23'].append(1.23-dG23); ORR['dg34'].append(1.23-dG12); ORR['dg41'].append(1.23-dG41)
    ORR['overP'].append(overP_orr); ORR['onsetP'].append(onsetP_orr)
    
def overpotential_oer(int1, int2, int3, int4, df, OER):
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
    int1, int2, int3, int4 = ints
    
    dG12 = df.loc[int2, 'dG'] - df.loc[int1, 'dG']
    dG23 = df.loc[int3, 'dG'] - df.loc[int2, 'dG']
    dG34 = df.loc[int4, 'dG'] - df.loc[int3, 'dG']
    dG41 = 4.92 - dG12 - dG23 - dG34
    if any(np.isnan(value) for value in [dG12, dG23, dG34, dG41]):
        onsetP = np.nan
        overP = np.nan
    else:
        onsetP = max(dG12, dG23, dG34, dG41)
        overP = onsetP - 1.23
    OER['int1'].append(int1); OER['int2'].append(int2); OER['int3'].append(int3); OER['int4'].append(int4)
    OER['dg12'].append(dG12); OER['dg23'].append(dG23); OER['dg34'].append(dG34); OER['dg41'].append(dG41)
    OER['overP'].append(overP); OER['onsetP'].append(onsetP)
    
def overpotential_orr(int1, int2, int3, int4, df, ORR):
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
    int1, int2, int3, int4 = ints
    
    dG12 = df.loc[int2, 'dG'] - df.loc[int1, 'dG']
    dG23 = df.loc[int3, 'dG'] - df.loc[int2, 'dG']
    dG34 = df.loc[int4, 'dG'] - df.loc[int3, 'dG']
    dG41 = -4.92 - dG12 - dG23 - dG34
    if any(np.isnan(value) for value in [dG43, dG32, dG34, dG41]):
        onsetP = np.nan
        overP = np.nan
    else:
        onsetP = -max(dG12, dG23, dG34, dG41)
        overP = 1.23 + max(dG12, dG23, dG34, dG41)
    ORR['int1'].append(int1); ORR['int2'].append(int2); ORR['int3'].append(int3); ORR['int4'].append(int4)
    ORR['dg12'].append(dG12); ORR['dg23'].append(dG23); ORR['dg34'].append(dG34); ORR['dg41'].append(dG41)
    ORR['overP'].append(overP); ORR['onsetP'].append(onsetP)
    
for dir in dirs:
    os.chdir(dir)
    print(dir)
    basename = os.path.basename(os.path.normpath(dir))
    A, B = basename.split('_', 1)
    df = pd.DataFrame()
    OER = {'int1': [], 'int2': [], 'int3': [], 'int4': [], 'dg12': [], 'dg23': [], 'dg34': [], 'dg41': [], 'overP': [], 'onsetP': []}
    ORR = {'int1': [], 'int2': [], 'int3': [], 'int4': [], 'dg12': [], 'dg23': [], 'dg34': [], 'dg41': [], 'overP': [], 'onsetP': []}
    
    json_path = '/pscratch/sd/j/jiuy97/6_MNC/empty/2_/final_with_calculator.json'
    atoms = read(json_path)
    df.loc['vac', 'E'] = atoms.get_potential_energy()
    
    # Iterate through each main directory to extract E0 values and plot
    for main_dir in main_dirs:
        min_e0 = get_energy(main_dir, sub_dirs)
    
        if min_e0 is None:
            # print(f"Missing data in directory '{main_dir}' for plotting.")
            df.loc[main_dir, 'E'] = np.nan
            continue
        else:
            df.loc[main_dir, 'E'] = min_e0
                
    charge = elements_data[B]['cation_charge']
    cation = metal_df.at[B, 'energy'] + charge * elements_data[B]['electrode_potential']
    
    df.loc['vac', ['#H', '#O', '#OH', '#OOH']] = [2, 0, 0, 0]
    df.loc['clean', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 0] # [energy, #Hs, #Os, #OHs, #OOHs]
    df.loc['mh', ['#H', '#O', '#OH', '#OOH']] = [1, 0, 0, 0]
    df.loc['nh', ['#H', '#O', '#OH', '#OOH']] = [1, 0, 0, 0]
    df.loc['o', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 0]
    df.loc['oh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 0]
    df.loc['ohoh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 2, 0]
    df.loc['oh-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 2, 0]
    df.loc['ohooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['oohoh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['oh-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['ooh-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 1]
    df.loc['oho', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['oh-o', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['o-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['oo', ['#H', '#O', '#OH', '#OOH']] = [0, 2, 0, 0]
    df.loc['o-o', ['#H', '#O', '#OH', '#OOH']] = [0, 2, 0, 0]
    df.loc['oooh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['ooho', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['o-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['ooh-o', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['oohooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 2]
    df.loc['ooh-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 2]

    df['G'] = df['E'] + dgh * df['#H'] + dgo * df['#O'] + dgoh * df['#OH'] + dgooh * df['#OOH']
    df['dG'] = df['G'] - df.loc['clean', 'E'] - gh * df['#H'] - go * df['#O'] - goh * df['#OH'] - gooh * df['#OOH']
    df.loc['vac', 'dG'] += cation
    
    overpotential('clean', 'oh', 'o', 'ooh', df, OER, ORR)
    if A == '1' and B == 'Fe':
        overpotential('oh', 'oh-oh', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, OER, ORR)
        overpotential('ooh', ('oh-ooh', 'ooh-oh'), ('ooh-o', 'o-ooh'), 'ooh-ooh', df, OER, ORR)
        overpotential('oh', 'o', 'ooh', ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('oh', 'o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', 'ooh', ('ooh-oh', 'oh-ooh'), ('ooh-o', 'o-ooh'), df, OER, ORR)        
        overpotential('o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), ('ooh-o', 'o-ooh'), df, OER, ORR)
    elif A == '2' and B == 'Co':
        overpotential('oh', 'oh-oh', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, OER, ORR)
        overpotential('ooh', ('oh-ooh', 'ooh-oh'), ('ooh-o', 'o-ooh'), 'ooh-ooh', df, OER, ORR)
        overpotential('oh', 'o', 'ooh', ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('oh', 'o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, OER, ORR)
        overpotential('o', 'ooh', ('ooh-oh', 'oh-ooh'), ('ooh-o', 'o-ooh'), df, OER, ORR)        
        overpotential('o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), ('ooh-o', 'o-ooh'), df, OER, ORR)
    elif A == '3' and B == 'Mo':
        overpotential('oh', 'o', 'oho', ('oohoh', 'ohooh'), df, OER, ORR)
        overpotential('o', 'oho', 'oo', ('oooh', 'ooho'), df, OER, ORR)
        overpotential('oh', 'ohoh', 'oho', 'oo', df, OER, ORR)
        overpotential('oh', 'o', 'oho', 'oo', df, OER, ORR)
        
    # Define surfaces with extracted E0 values
    surfs = [
        df.loc['vac', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        df.loc['o-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
    ]
    
    if A == '3' and B == 'Mo':
        surfs.append(df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    # if A == '1' and B == 'Fe':
    #     # surfs.append(df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     # surfs.append(df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    # elif A == '2' and B == 'Co':
    #     # surfs.append(df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     # surfs.append(df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    # elif A == '3' and B == 'Mo':
    #     surfs.append(df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    #     surfs.append(df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist())

    
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
    
    # print(lowest_surfaces)
    # print(U2)
    # print(len(lowest_surfaces))
    # print(len(U2))
    # print(surfs)
    # print(uniquesurf)
    # print(crossover)
    
    # plt.clf()
    # fig = plt.figure(figsize=fig_size, dpi=300)
    # ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    # ax.axis([0, 14, Umin, Umax])
    # ax.set_xlabel(r'pH', fontsize='large')
    # ax.set_ylabel(r'U/V', fontsize='large')
    # current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    # extraticks = [1.23]
    # combined_ticks = sorted(set(current_yticks) | set(extraticks))
    # plt.yticks(combined_ticks)
    # plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # for i in range(len(uniquesurf)):
    #     k = uniquesurf[i]
    #     label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
    #     plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
    #                      facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
    #     plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)
    # plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    # if A == '1' and B == 'Fe':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][1] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.70, 
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][1] - 0.96, # 1.14, 
    #             r"S$_4$$\rightarrow$S$_7$$\rightarrow$S$_{10}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][7]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    # elif A == '2' and B == 'Co':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][3] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][3] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.94, 
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][3] - 0.74, 
    #             r"S$_8$$\rightarrow$S$_{12}$$\rightarrow$S$_{13}$$\rightarrow$S$_{14}$: " + f"{OER['overP'][3]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    # elif A == '3' and B == 'Mo':
    #     plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.8, OER['onsetP'][2] - 1.30, 
    #             r"S$_5$$\rightarrow$S$_8$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$: " + f"{OER['overP'][2]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
    #            ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
    #            fancybox=True, shadow=True)
    # plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_pourbaix_oer.png', bbox_inches='tight')
    # print(f"Figure saved as {A}{B}_pourbaix_oer.png")
    # plt.close()

    # plt.clf()
    # fig = plt.figure(figsize=fig_size, dpi=300)
    # ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    # ax.axis([0, 14, Umin, Umax])
    # ax.set_xlabel(r'pH', fontsize='large')
    # ax.set_ylabel(r'U/V', fontsize='large')
    # current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    # extraticks = [1.23]
    # combined_ticks = sorted(set(current_yticks) | set(extraticks))
    # plt.yticks(combined_ticks)
    # plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # for i in range(len(uniquesurf)):
    #     k = uniquesurf[i]
    #     label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
    #     plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
    #                      facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
    #     plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)
    # plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    # if A == '1' and B == 'Fe':
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_7$$\rightarrow$S$_4$$\rightarrow$S$_3$$\rightarrow$S$_0$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    # elif A == '2' and B == 'Co':
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     # plt.plot(pH2, ORR['onsetP'][5] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_7$$\rightarrow$S$_4$$\rightarrow$S$_3$$\rightarrow$S$_0$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    # elif A == '3' and B == 'Mo':
    #     plt.plot(pH2, ORR['onsetP'][2] - pH2 * const, '--', color='brown', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][2] - 0.34,
    #             r"S$_{10}$$\rightarrow$S$_9$$\rightarrow$S$_7$$\rightarrow$S$_4$: " + f"{ORR['overP'][2]:.2f} eV", 
    #             color='brown', rotation=-9.5, fontsize=10)
    # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
    #            ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
    #            fancybox=True, shadow=True)
    # plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_pourbaix_orr.png', bbox_inches='tight')
    # print(f"Figure saved as {A}{B}_pourbaix_orr.png")
    # plt.close()
    
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
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                         facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)
    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    if A == '1' and B == 'Fe':
        plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
        plt.plot(pH2, OER['onsetP'][1] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][0] - 0.70, 
                r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
                color='black', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][1] - 0.96, # 1.14, 
                r"S$_4$$\rightarrow$S$_7$$\rightarrow$S$_{10}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][1]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][0] - 0.58, 
                r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
                color='gray', rotation=-9.5, fontsize=10)
    elif A == '2' and B == 'Co':
        plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
        plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][0] - 0.70,
                r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
                color='black', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][2] - 0.98, 
                r"S$_5$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][2]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][0] - 0.58, 
                r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
                color='gray', rotation=-9.5, fontsize=10)
    elif A == '3' and B == 'Mo':
        plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][2] - pH2 * const, '--', color='green', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.8, OER['onsetP'][2] - 1.30, 
                r"S$_5$$\rightarrow$S$_8$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$: " + f"{OER['overP'][2]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][2] - 0.34,
                r"S$_{11}$$\rightarrow$S$_{10}$$\rightarrow$S$_8$$\rightarrow$S$_5$: " + f"{ORR['overP'][2]:.2f} eV", 
                color='green', rotation=-9.5, fontsize=10)
    plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_pourbaix_full.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_full.png")
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
    if A == '1' and B == 'Fe':
        plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
        plt.plot(pH2, OER['onsetP'][1] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][0] - 0.70, 
                r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
                color='black', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][1] - 0.96, # 1.14, 
                r"S$_4$$\rightarrow$S$_7$$\rightarrow$S$_{10}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][1]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][0] - 0.58, 
                r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
                color='gray', rotation=-9.5, fontsize=10)
    elif A == '2' and B == 'Co':
        plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
        plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][0] - 0.70,
                r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
                color='black', rotation=-9.5, fontsize=10)
        ax.text(6.5, OER['onsetP'][2] - 0.98, 
                r"S$_5$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][2]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][0] - 0.58, 
                r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
                color='gray', rotation=-9.5, fontsize=10)
    elif A == '3' and B == 'Mo':
        plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
        plt.plot(pH2, ORR['onsetP'][2] - pH2 * const, '--', color='green', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(6.8, OER['onsetP'][2] - 1.30, 
                r"S$_5$$\rightarrow$S$_8$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$: " + f"{OER['overP'][2]:.2f} eV", 
                color='red', rotation=-9.5, fontsize=10)
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(0.2, ORR['onsetP'][2] - 0.34,
                r"S$_{11}$$\rightarrow$S$_{10}$$\rightarrow$S$_{8}$$\rightarrow$S$_{5}$: " + f"{ORR['overP'][2]:.2f} eV", 
                color='green', rotation=-9.5, fontsize=10)
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_pourbaix_clean.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_clean.png")
    plt.close()
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([-1.0, 2.5, -600, 200])
    # if A=='3' and B=='Mo':
    #     ax.axis([-1.0, 2.5, -900, 300])
    # else:
    #     ax.axis([-1.0, 2.5, -600, 200])
    ax.set_xlabel(r'RHE (V)', fontsize='large')
    ax.set_ylabel(r'$\Delta$G (kJ/mol)', fontsize='large')
    xx = np.arange(-1.00, 2.55, 0.01)
    for k in range(nsurfs):
        label = r"S$_{%i}$(H: %i O: %i OH: %i OOH: %i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        dg_value = dg(k, 0, xx)
        if dg_value is not None:
            ax.plot(xx, dg_value * kjmol, '-', lw=1, c=color[k], label=label)
        else:
            print(f"Skipping plot for surface {k} due to missing data.")
    plt.xlim(-1.0, 2.5)
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
    #            ncol=2, columnspacing=1.0, labelspacing=0.3, handlelength=2, fontsize=10,
    #            fancybox=True, shadow=True)
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_pourbaix.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix.png")
    # plt.show()
    plt.close()

    df['#H'] = df['#H'].astype(int)
    df['#O'] = df['#O'].astype(int)
    df['#OH'] = df['#OH'].astype(int)
    df['#OOH'] = df['#OOH'].astype(int)
    df['E'] = df['E'].round(2)
    df['G'] = df['G'].round(2)
    df['dG'] = df['dG'].round(2)
    df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_energies.tsv', sep='\t') #, index=False)
    print(f"Data saved as {A}{B}_energies.png")

    OER_df = pd.DataFrame(OER)
    OER_df['dg12'] = OER_df['dg12'].round(2)
    OER_df['dg23'] = OER_df['dg23'].round(2)
    OER_df['dg34'] = OER_df['dg34'].round(2)
    OER_df['dg41'] = OER_df['dg41'].round(2)
    OER_df['overP'] = OER_df['overP'].round(2)
    OER_df['onsetP'] = OER_df['onsetP'].round(2)
    OER_df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_oer.tsv', sep='\t') #, index=False)
    print(f"Data saved as {A}{B}_oer.png")

    ORR_df = pd.DataFrame(ORR)
    ORR_df['dg12'] = ORR_df['dg12'].round(2)
    ORR_df['dg23'] = ORR_df['dg23'].round(2)
    ORR_df['dg34'] = ORR_df['dg34'].round(2)
    ORR_df['dg41'] = ORR_df['dg41'].round(2)
    ORR_df['overP'] = ORR_df['overP'].round(2)
    ORR_df['onsetP'] = ORR_df['onsetP'].round(2)
    ORR_df.to_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/{A}{B}_orr.tsv', sep='\t') #, index=False)
    print(f"Data saved as {A}{B}_orr.png")
    