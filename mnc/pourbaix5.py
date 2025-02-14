import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
   
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

color = [
    # 'cornflowerblue', ## Sv
    'lightsteelblue', ## S0
    'darkgray', ## S1
    'yellowgreen', # S2
    'red', # S3
    'tan', ## S4
    'pink', ## S5
    'green', # S6
    'salmon', ## S7  
    'blue', # S8
    'khaki', ## S9
    'salmon', ## S10
    'plum', ## S11
]
pH2 = np.arange(0, 14.01, 0.01)

# path = '/Users/hailey/Desktop/jan20/pourbaix/'
# metal_path = '/Users/hailey/Desktop/jan20/metals.tsv'
path = '/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/'
metal_path = '/pscratch/sd/j/jiuy97/6_MNC/figures/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
metals = ['Fe', 'Co', 'Mo']

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
                + bulk_metal
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
  
for i, metal in enumerate(metals):
    
    A, B = i+1, metal                
    bulk_metal = metal_df.at[B, 'energy']
    
    df = pd.DataFrame()
    df_path = os.path.join(path, f'{A}{B}_energies.tsv')
    df = pd.read_csv(df_path, delimiter='\t', index_col=0)
    df.loc['vac', 'E'] = -.28183609E+03            
    df.loc['vac', ['#H', '#O', '#OH', '#OOH']] = [2, 0, 0, 0]
    
    OER = pd.DataFrame()
    OER_path = os.path.join(path, f'{A}{B}_oer.tsv')
    OER = pd.read_csv(OER_path, delimiter='\t', index_col=0)

    ORR = pd.DataFrame()
    ORR_path = os.path.join(path, f'{A}{B}_orr.tsv')
    ORR = pd.read_csv(ORR_path, delimiter='\t', index_col=0)
    
    surfs = [
        df.loc['vac', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
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
    
    surfs = [surf for surf in surfs if not any(pd.isna(x) for x in surf)]
    
    nsurfs = len(surfs)
    lowest_surfaces = []
    
    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs) if dg(k, 0, j) is not None]
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
    ax.axis([-1.0, 2.5, -600, 200])
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
    plt.savefig(f'{path}{A}{B}_pourbaix.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix.png")
    # plt.show()
    plt.close()
    