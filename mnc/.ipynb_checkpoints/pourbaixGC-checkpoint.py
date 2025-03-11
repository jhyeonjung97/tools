import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
calmol = 23.061
pHrange = np.arange(0, 14.1, 0.05)
Umin, Umax = -1.0, 3.0
Urange = np.arange(Umin, Umax + 0.06 * 14, 0.05)

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

water = 2.4583 # the standard Gibbs free energy of formation of water

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
        most_stable_dir = f"{root}/pourbaix/1_Fe/{main_dir}/most_stable"
    done_path = os.path.join(most_stable_dir, "DONE")
    json_path = os.path.join(most_stable_dir, "final_with_calculator.json")
    if os.path.isfile(done_path) and os.path.isfile(json_path):
        atoms = read(json_path)
        return atoms.get_potential_energy()
    return None

def dg_surf(k, pH, U):
    if surfs[k][0] is None:
        return None
    dg = (
        (surfs[k][6]*(U**2) + surfs[k][7]*U + surfs[k][8])
        - (surfs[0][6]*(U**2) + surfs[0][7]*U + surfs[0][8])
        + surfs[k][2] * (dgh -gh + 1 * (U + pH * const))
        + surfs[k][3] * (dgo -go - 2 * (U + pH * const))
        + surfs[k][4] * (dgoh -goh - (U + pH * const))
        + surfs[k][5] * (dgooh -gooh - 3 * (U + pH * const))
    )
    if all(x != 0 for x in [surfs[k][6], surfs[k][7], surfs[k][8], surfs[0][6], surfs[0][7], surfs[0][8]]):
        print(surfs[k][6], surfs[k][7], surfs[k][8], surfs[0][6], surfs[0][7], surfs[0][8])
    return dg

def dg_ion(k, pH, U):
    dg = (
        surfs[k][0]
        - (surfs[0][6]*(U**2) + surfs[0][7]*U + surfs[0][8])
        + 2 * (dgh -gh + 1 * (U + pH * const))
        + surfs[k][2] * (1 * (U + pH * const))
        + surfs[k][3] * (-2 * (U + pH * const))
        + surfs[k][4] * (-1 * (U + pH * const))
        + surfs[k][5] * (-3 * (U + pH * const))
        - surfs[k][1] * U
    )
    return dg

df = pd.DataFrame()
json_path = f'{root}/empty/2_/final_with_calculator.json'
atoms = read(json_path)
vac = atoms.get_potential_energy()

df.loc['Fe', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
df.loc['Fe²⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-20.300/calmol, +2, 0, 0, 0, 0, 0, 0, 0]
df.loc['HFeO²⁻', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-90.627/calmol, -2, 1, 1, 0, 0, 0, 0, 0]
df.loc['Fe³⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-2.530/calmol, +3, 0, 0, 0, 0, 0, 0, 0]
df.loc['FeOH²⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-55.910/calmol, +2, 1, 1, 0, 0, 0, 0, 0]
df.loc['Fe(OH)₂⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-106.200/calmol, +1, 2, 2, 0, 0, 0, 0, 0]
df.loc['FeO', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-58.880/calmol, 0, 0, 1, 0, 0, 0, 0, 0]
# df.loc['Fe₃O₄', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-242.400/calmol, 0, 0, 0, 0, 0, 0, 0, 0]
# df.loc['Fe₂O₃', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-177.100/calmol, 0, 0, 0, 0, 0, 0, 0, 0]
df.loc['Fe(OH)₂', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-115.570/calmol, 0, 2, 2, 0, 0, 0, 0, 0]
df.loc['Fe(OH)₃', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [-166.000/calmol, 0, 3, 3, 0, 0, 0, 0, 0]
        
df['E'] += vac + bulk_metal + water * df['#O']
df['#H'] += 2

df.loc['vac', 'E'] = vac
for main_dir in main_dirs:
    min_e0 = get_energy(main_dir)
    if min_e0 is None:
        df.loc[main_dir, 'E'] = np.nan
    else:
        df.loc[main_dir, 'E'] = min_e0

df.loc['vac', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 2, 0, 0, 0, -0.3342, -0.1079, 0]
df.loc['clean', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, 0, -0.3660, 0.0665, 0]
df.loc['mh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 0, 0, -0.3714, 0.2173, 0]
df.loc['nh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 1, 0, 0, 0, -0.3438, -0.3362, 0]
df.loc['o', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, 0, -0.3608, -0.6948, 0]
df.loc['oh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, 0, -0.8300, -0.4206, 0]
df.loc['ohoh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 2, 0, -0.7036, 0.3162, 0]
df.loc['oh-oh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 2, 0, -0.7214, 0.3152, 0]
df.loc['ohooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, 1, -0.9681, -0.9980, 0]
df.loc['oohoh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, 1, -0.9971, -1.1225, 0]
df.loc['oh-ooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, 1, -0.8134, 0.0605, 0]
df.loc['ooh-oh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 1, 1, -0.9297, -0.7491, 0]
df.loc['ooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, 1, -0.8300, -0.4206, 0]
df.loc['oho', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, 0, -0.8385, -0.6094, 0]
df.loc['oh-o', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, 0, -0.7100, 0.1638, 0]
df.loc['o-oh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 1, 0, -0.7247, 1.7570, 0]
df.loc['oo', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 2, 0, 0, -0.6536, 2.1095, 0]
df.loc['o-o', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 2, 0, 0, -0.6429, 2.1679, 0]
df.loc['oooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, 1, -0.6918, 1.3097, 0]
df.loc['ooho', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, 1, -0.7158, 0.8580, 0]
df.loc['o-ooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, 1, -0.7384, 2.0032, 0]
df.loc['ooh-o', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 1, 0, 1, -0.6788, 1.2123, 0]
df.loc['oohooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, 2, -0.6258, 1.0597, 0]
df.loc['ooh-ooh', ['#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']] = [0, 0, 0, 0, 2, -0.5632, 1.0037, 0]
df.at['vac', 'E'] += bulk_metal
df['C'] = df['E']
df[['A', 'B']] = 0

print(df)

surfs = [
    df.loc['vac', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['clean', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['mh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['nh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['oh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['oh-oh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['oh-o', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    df.loc['o-oh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['o-o', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['o', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    df.loc['ohoh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['ohooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oohoh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oh-ooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['ooh-oh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oho', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['oo', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['oooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooho', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['o-ooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh-o', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['oohooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    # df.loc['ooh-ooh', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),  
    
    # df.loc['Fe²⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['HFeO²⁻', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #1
    # df.loc['Fe³⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #2
    # df.loc['FeOH²⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #3
    # df.loc['Fe(OH)₂⁺', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(), #4
    # df.loc['FeO', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # # df.loc['Fe₃O₄', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # # df.loc['Fe₂O₃', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['Fe(OH)₂', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
    # df.loc['Fe(OH)₃', ['E', '#e', '#H', '#O', '#OH', '#OOH', 'A', 'B', 'C']].tolist(),
]
surfs = [surf for surf in surfs if not any(pd.isna(x) for x in surf)]
nsurfs = len(surfs)

lowest_surfaces = np.ones((len(Urange),len(pHrange)))*-1
    
pHindex = 0
for pH in pHrange:
    Uindex = 0
    for U in Urange:
        values = []
        for k, surf in enumerate(surfs):
            if surf[1] != 0:
                values.append(dg_ion(k, pH, U))
            else:
                values.append(dg_surf(k, pH, U))
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Uindex][pHindex] = sorted_values[0]
        Uindex+=1
    pHindex+=1
    
fig = plt.figure(figsize=(8, 6), dpi=100)
ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])

# Set Axes Limits and Labels
ax.axis([0, 14, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

pH, U = np.meshgrid(pHrange, Urange)
colors = plt.cm.tab20.colors
cmap = mcolors.ListedColormap(colors)

for k in range(nsurfs):
    if k in lowest_surfaces:
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][2], surfs[k][3], surfs[k][4], surfs[k][5])
        plt.plot([], [], color=colors[k], linewidth=5, label=label)
        
unique_surfs = np.unique(lowest_surfaces)
selected_colors = [colors[int(k) % len(colors)] for k in unique_surfs]  
lowest_cmap = mcolors.ListedColormap(selected_colors)

selected_surfaces = []
for i in lowest_surfaces.flatten(order='F'):
    selected_surfaces.append(i)
unique_surfaces = np.unique(selected_surfaces)
selected_colors = [colors[int(k) % len(colors)] for k in unique_surfaces]  
lowest_cmap = mcolors.ListedColormap(selected_colors)

print(surfs)
print(lowest_surfaces)
# print(selected_surfaces)
print(unique_surfaces)

plt.pcolormesh(pH, U, lowest_surfaces, shading='auto', cmap=lowest_cmap, vmin=0, vmax=nsurfs-1)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1,
       fontsize='x-small', handlelength=3, edgecolor='black')

plt.plot(pHrange, 1.23-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 1.00, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$',
        color='blue', rotation=-8, fontsize=10)
plt.plot(pHrange, 0-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, -0.15 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$',
        color='blue', rotation=-8, fontsize=10)

fig.savefig(f'{root}/figures/pourbaixGC.png', bbox_inches='tight') 