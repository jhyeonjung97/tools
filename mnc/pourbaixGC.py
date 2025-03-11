import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

root = '/pscratch/sd/j/jiuy97/6_MNC'
main_dirs = ["clean", "mh", "nh", "oh", "o",
             "ohoh", "oh-oh", "ohooh", "oohoh", "oh-ooh", "ooh-oh",
             "ooh", "oho", "oh-o", "o-oh", "oo", "o-o",
             "oooh", "ooho", "o-ooh", "ooh-o", "oohooh", "ooh-ooh",]
             # "oo", "oo-oo", "oo-ohh", "oo-o", "oo-oh", "oo-ooh",
             # "ohh", "ohh-oo", "ohh-ohh", "ohh-o", "ohh-oh", "ohh-ooh",]

# constants
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485
pHrange = np.arange(0, 14.01, 0.01)
Umin, Umax = -0.5, 2.5 + 0.06 * 14
Urange = np.arange(Umin, Umax, 0.01)

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

figure_path = f'{root}/figures/pourbaix'
metal_path = '~/bin/tools/mnc/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
bulk_metal = metal_df.at['Fe', 'energy']

def get_energy(main_dir):
    if main_dir == 'clean':
        most_stable_dir = f"{root}/0_clean/3d/6_Fe/most_stable/relaxed"
    elif main_dir == 'o':
        most_stable_dir = f"{root}/1_O/2_Fe/most_stable/relaxed"
    elif main_dir == 'oh':
        most_stable_dir = f"{root}/2_OH/2_Fe/most_stable/relaxed"
    elif main_dir == 'ooh':
        most_stable_dir = f"{root}/3_OOH/2_Fe/most_stable/relaxed"
    else:
        most_stable_dir = os.path.join(main_dir, "most_stable")
    done_path = os.path.join(most_stable_dir, "DONE")
    json_path = os.path.join(most_stable_dir, "final_with_calculator.json")
    if os.path.isfile(done_path) and os.path.isfile(json_path):
        atoms = read(json_path)
        return atoms.get_potential_energy()
    return None
    
def addH(pH, U):
    return -gh + dgh + 1 * (U + pH * const)
    
def addO(pH, U):
    return -go + dgo - 2 * (U + pH * const)

def addOH(pH, U):
    return -goh + dgoh - (U + pH * const)

def addOOH(pH, U):
    return -gooh + dgooh - 3 * (U + pH * const)

def dg(i, pH, U):
    if surfs[i][0] is None:
        return None
    dg = (
        (surfs[i][5]*(U**2) + surfs[i][6]*U + surfs[i][7])
        - (surfs[0][5]*(U**2) + surfs[0][6]*U + surfs[0][7])
        + surfs[i][1] * addH(pH, U) 
        + surfs[i][2] * addO(pH, U) 
        + surfs[i][3] * addOH(pH, U) 
        + surfs[i][4] * addOOH(pH, U)
    )
    if i == 0 and surfs[i][1] == 2:
        return dg + bulk_metal
    return dg
    
df = pd.DataFrame()
json_path = f'{root}/empty/2_/final_with_calculator.json'
atoms = read(json_path)
df.loc['vac', 'E'] = atoms.get_potential_energy()

for main_dir in main_dirs:
    min_e0 = get_energy(main_dir)
    if min_e0 is None:
        df.loc[main_dir, 'E'] = np.nan
    else:
        df.loc[main_dir, 'E'] = min_e0

df.loc['vac', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [2, 0, 0, 0, -0.3342, -0.1079, -421.8599]
df.loc['clean', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, -0.3660, 0.0665, -441.1985]
df.loc['mh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [1, 0, 0, 0, -0.3714, 0.2173, -429.6398]
df.loc['nh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [1, 0, 0, 0, -0.3438, -0.3362, -452.2326]
df.loc['o', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 0, -0.3608, -0.6948, -431.9379]
df.loc['oh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, -0.8300, -0.4206, -592.8774]
df.loc['ohoh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 2, 0, -0.7036, 0.3162, -633.5199]
df.loc['oh-oh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 2, 0, -0.7214, 0.3152, -610.2184]
df.loc['ohooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, -0.9681, -0.9980, -607.6129]
df.loc['oohoh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, -0.9971, -1.1225, -612.2014]
df.loc['oh-ooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, -0.8134, 0.0605, -632.3953]
df.loc['ooh-oh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, -0.9297, -0.7491, -653.2411]
df.loc['ooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, -0.8300, -0.4206, -592.8774]
df.loc['oho', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 1, 0, -0.8385, -0.6094, -648.7112]
df.loc['oh-o', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 1, 0, -0.7100, 0.1638, -627.2593]
df.loc['o-oh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 1, 0, -0.7247, 1.7570, -441.1449]
df.loc['oo', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 2, 0, 0, -0.6536, 2.1095, -461.8944]
df.loc['o-o', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 2, 0, 0, -0.6429, 2.1679, -450.4421]
df.loc['oooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 1, -0.6918, 1.3097, -470.9077]
df.loc['ooho', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 1, -0.7158, 0.8580, -450.2448]
df.loc['o-ooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 1, -0.7384, 2.0032, -461.1960]
df.loc['ooh-o', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 1, -0.6788, 1.2123, -445.9051]
df.loc['oohooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 2, -0.6258, 1.0597, -445.8280]
df.loc['ooh-ooh', ['#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 2, -0.5632, 1.0037, -445.8658]

surfs = [
    df.loc['vac', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #1
    df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #2
    df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #3
    df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #4
    df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #5
    df.loc['o-oh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #6
    df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #7
    df.loc['o', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #8
    df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #9
    # df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #10
    # df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
]
surfs = [surf for surf in surfs if not any(pd.isna(x) for x in surf)]
nsurfs = len(surfs)
lowest_surfaces = []

# Loop find the coverage with the lowest dG at each value of pH, voltage
pHindex = 0
for pH in pHrange:
    Uindex = 0
    for U in Urange:
        values = []
        for k in range(nsurfs):
            values.append(dg(k, pH, U))
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Uindex][pHindex] = sorted_values[0]
        Uindex+=1
    pHindex+=1
    
fig = plt.figure(figsize=(8, 6), dpi=100)
ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])

cmapName = 'RdYlBu'
pH, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH, U, lowest_surfaces, shading='auto', norm=None, cmap=cmapName, alpha=0.85, vmin=0, vmax=nsurfs)
cmap = plt.get_cmap(cmapName, nsurfs+1)

for k in range(nsurfs): 
    label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
    plt.plot([], [], color=cmap(k), linewidth=5, label=label)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1,
       fontsize='x-small', handlelength=3, edgecolor='black')

plt.plot(pHrange, 1.23-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 0.95, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$',
        color='blue', rotation=-14, fontsize=10)
plt.plot(pHrange, 0-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, -0.2 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$',
        color='blue', rotation=-14, fontsize=10)

fig.savefig(f'{root}/figures/pourbaixGC.png', bbox_inches='tight') 