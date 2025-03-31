### Hailey ###
### Date: 2025-03-13 ###

import os
import numpy as np
from math import log10
from matplotlib import colormaps
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

GCDFT = False
# bulk_metal = -5.627328705 # Ti, eV
bulk_metal = -5.7426 # Ti, eV

# constants
kb = 8.617e-5 # eV/K
T = 298.15 # K
const = kb * T * np.log(10) # 0.0592 eV
water = 2.4583 # the standard Gibbs free energy of formation of water

# units
kjmol = 96.485
calmol = 23.061

# ticks
tick = 0.01
pHrange = np.arange(0, 14.1, tick)
Umin, Umax = -1.0, 3.0
Urange = np.arange(Umin, Umax + 0.06 * 14, tick)

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

def dg(k, pH, U, concentration):
    if GCDFT:
        surface_term = ((surfs[k][8]*(U**2) + surfs[k][9]*U + surfs[k][10])
                          - (surfs[0][8]*(U**2) + surfs[0][9]*U + surfs[0][10]))
    else:
        surface_term = surfs[k][0] - surfs[0][0]
    U_coeff = 1*surfs[k][4] - 1*surfs[k][5] - 2*surfs[k][6] - 3*surfs[k][7] - surfs[k][3]
    pH_coeff = 1*surfs[k][4] - 1*surfs[k][5] - 2*surfs[k][6] - 3*surfs[k][7]
    dg = surface_term + U_coeff * U + pH_coeff * const * pH
    if '(aq)' in surfs[k][11]:
        dg += const * log10(concentration)
    return dg
    
ions = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [ -75.100/calmol, 1, 0, +2, 0, 0, 0, 0, 0, 0, 0, 'Ti²⁺(aq)'],
    [ -83.600/calmol, 1, 0, +3, 0, 0, 0, 0, 0, 0, 0, 'Ti³⁺(aq)'],
    [-138.000/calmol, 1, 0, +2, 0, 0, 0, 0, 0, 0, 0, 'TiO²⁺(aq)'],
    [-228.460/calmol, 1, 0, -1, 1, 0, 3, 0, 0, 0, 0, 'HTiO₃⁻(aq)'],
    [-111.670/calmol, 1, 0, +2, 0, 0, 2, 0, 0, 0, 0, 'TiO₂²⁺(aq)'],
    [ -19.100/calmol, 0, 1, +0, 0, 0, 0, 0, 0, 0, 0, 'HNO₃(aq)'],
    [ -63.050/calmol, 0, 1, +0, 0, 0, 0, 0, 0, 0, 0, 'NHOH'],
    [ -19.000/calmol, 0, 1, +0, 0, 0, 0, 0, 0, 0, 0, 'NH+'],
    [  30.560/calmol, 0, 2, +0, 0, 0, 0, 0, 0, 0, 0, 'N2H4'],
    [  22.500/calmol, 0, 2, +0, 0, 0, 0, 0, 0, 0, 0, 'N2H4N2++'],
    [  -5.600/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NH2OH'],
    [ -13.540/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NH2OH2+'],
    [  71.300/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'HN3'],
    [  77.700/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N3-'],
    [   2.994/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N2'],
    [   8.600/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'H2N2O2'],
    [  33.200/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N2O2--'],
    [  18.200/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NH2O2-'],
    [ -12.820/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'HNO2'],
    [  -8.250/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NO2-'],
    [ -26.430/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'HNO3'],
    [ -26.430/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NO3-'],
]

solids = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0,               1, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'Ti(s)'],
    [-116.920/calmol, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, 'TiO'],
    [-342.310/calmol, 2, 0, +0, 0, 0, 3, 0, 0, 0, 0, 'Ti₂O₃'],
    [-250.905/calmol, 1, 0, +0, 0, 3, 0, 0, 0, 0, 0, 'Ti(OH)₃'],
    [-553.120/calmol, 3, 0, +0, 0, 0, 5, 0, 0, 0, 0, 'Ti₃O₅'],
    [-212.330/calmol, 1, 0, +0, 0, 0, 2, 0, 0, 0, 0, 'TiO₂'],
    [-252.990/calmol, 1, 0, +0, 1, 1, 2, 0, 0, 0, 0, 'TiO₂·H₂O'],
    [  32.000/calmol, 0, 2, +0, 0, 0, 5, 0, 0, 0, 0, 'N₂O₅'],
]

gases = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-3.976/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NH3'],
    [78.500/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'HN3'],
    [81.476/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N'],
    [     0/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N2'],
    [24.760/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N2O'],
    [20.719/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NO'],
    [12.390/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'NO2'],
    [23.491/calmol, 0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'N2O4'],
]

nions, nsolids, ngases = len(ions), len(solids), len(gases)
for i in range(nions):
    ions[i][0] += water * (ions[i][5] + ions[i][6] + 2*ions[i][7]) + bulk_metal * ions[i][1] + gas_nitrogen * ions[i][2]
    if ions[i][1] > 1:
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]
    elif ions[i][2] > 1:
        ions[i] = [x / ions[i][2] if isinstance(x, (int, float)) else x for x in ions[i]]
    ions[i][10] = ions[i][0]

for s in range(nsolids):
    solids[s][0] += water * (solids[s][5] + solids[s][6] + 2*solids[s][7]) + bulk_metal * solids[s][1] + gas_nitrogen * solids[s][2]
    if solids[s][1] > 1:
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]
    elif solids[s][2] > 1:
        solids[s] = [x / solids[s][2] if isinstance(x, (int, float)) else x for x in solids[s]]
    solids[s][10] = solids[s][0]

for g in range(ngases):
    gases[g][0] += water * (gases[g][5] + gases[g][6] + 2*gases[g][7]) + bulk_metal * gases[g][1] + gas_nitrogen * gases[g][2]
    if gases[g][1] > 1:
        gases[g] = [x / gases[g][1] if isinstance(x, (int, float)) else x for x in gases[g]]
    elif gases[g][2] > 1:
        gases[g] = [x / gases[g][2] if isinstance(x, (int, float)) else x for x in gases[g]]
    gases[g][10] = gases[g][0]

surfs = [
    # ['E', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-269.569746,   0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'vac'],
    [-279.340765,   0, 0, +0, 2, 0, 0, 0, 0, 0, 0, 'vac(H₂)'],
    [-277.23249688, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'clean'],
    [-304.97333512, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₃(mono)'],
    [-305.01556231, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₃(mono)'], ##
    [-305.14584901, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₃(bi)'],
    [-300.39230369, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*OH+*OH(anti)'],
    [-301.21590570, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*OH+*OH(adj)'],
    [-297.04551538, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*HNO'],
    [-310.30495026, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₂(mono)+*OH(adj)'],
    [-300.07883828, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*HNOH'],
    [-283.58401076, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*H+*H(anti)'],
    [-284.73355007, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*H+*H(adj)'],
    [-294.10445417, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*H+*NO(anti)'],
    [-281.52877857, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*H'],
    [-284.81503764, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*N'],
    [-290.82616987, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NH'],
    [-295.02930066, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NH₂'],
    [-306.36204655, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NH₂+*OH(adj)'],
    [-292.62371142, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO'],
    [-299.05807974, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₂'],
    [-307.08538159, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₃(bi)+*H(anti)'],
    [-296.12169107, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NOH'],
    [-303.48187042, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*NO+*OH(adj)'],
    [-299.63088262, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*NO₂(bi)'],
    [-291.16160779, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*ON'],
    [-301.40150027, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*ONH₂'],
    [-306.42159787, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*O+*NH₃(adj)'],
    [-303.06694618, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*HNO₂'],
    [-295.88614166, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*ONH'],
    [-286.71126084, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*O'],
    [-290.31253022, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*O+*O(anti)'],
]

nsurfs = len(surfs)
ref0 = surfs[0][0]
ref9 = surfs[0][9]
for k in range(nsurfs):
    formation_energy_corr = (
        - surfs[k][3] * (gh - dgh) # H
        - surfs[k][4] * (goh - dgoh) # OH
        - surfs[k][5] * (go - dgo) # O
        - surfs[k][6] * (gooh - dgooh) # OOH
    )
    surfs[k][0] = surfs[k][0] - ref0 + formation_energy_corr 
    surfs[k][9] = surfs[k][9] - ref9 + formation_energy_corr 

new_surfs = []
for k in range(nsurfs):
    if surfs[k][1] == 0:
        for i in range(nions):
            new_surf = []
            for j in range(10):
                new_surf.append(surfs[k][j] + ions[i][j])
            new_surf.append(surfs[k][10] + '+' + ions[i][10])
            new_surfs.append(new_surf)
        for s in range(nsolids):
            new_surf = []
            for j in range(10):
                new_surf.append(surfs[k][j] + solids[s][j])
            new_surf.append(surfs[k][10] + '+' + solids[s][10])
            new_surfs.append(new_surf)

surfs.extend(new_surfs)  # Add new surfaces after looping
nsurfs = len(surfs)  # Update length

lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

for i in range(nsurfs):
    print(surfs[i][10], surfs[i][0])

pHindex = 0
for pH in pHrange:
    Uindex = 0
    for U in Urange:
        values = []
        for k in range(nsurfs):
            value = dg(k, pH, U, concentration=1e-6)
            values.append(value)
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Uindex][pHindex] = sorted_values[0]
        Uindex+=1
    pHindex+=1

min_coords = {}
n_rows, n_cols = lowest_surfaces.shape

for j in range(n_cols):   # loop over pH (columns)
    for i in range(n_rows):       # loop over U (rows)
        sid = int(lowest_surfaces[i, j])
        x = pHrange[j]   # pH varies along columns (x-axis)
        y = Urange[i]    # U varies along rows    (y-axis)
        if sid not in min_coords:
            min_coords[sid] = (x, y)
        else:
            current_x, current_y = min_coords[sid]
            if x < current_x or (x == current_x and y < current_y):
                min_coords[sid] = (x, y)

for sid in sorted(min_coords):
    x, y = min_coords[sid]
    print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}")
    
# Set Axes Limits and Labels
fig, ax = plt.subplots(figsize=(7, 6), dpi=100)
ax.axis([0, 14, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

# Check unique values in lowest_surfaces and create a list of unique surface IDs
unique_ids = np.unique(lowest_surfaces)
nsurfs = len(unique_ids)

# Define base color groups and their initial indices
base_colors = {
    'dodgerblue': 1,
    'orange': 1,
    'limegreen': 1,
    'gold': 1,
    'mediumblue': 1,
    'darkgoldenrod': 1,
    'hotpink': 1,
    'silver': 1,
}

# Generate custom colormaps and shades from white → base color
cmaps = {}
shades = {}

for base_color in base_colors:
    cmap = LinearSegmentedColormap.from_list(f"custom_{base_color}", ["white", base_color])
    cmaps[base_color] = cmap
    if base_color == 'dodgerblue' or base_color == 'orange':
        shades[base_color] = [cmap(i) for i in np.linspace(0, 1, 4)]
    else:
        shades[base_color] = [cmap(i) for i in np.linspace(0, 1, 4)]
        
# Map surface ID to corresponding color shade
color_mapping = {}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][10]
    if 'vac(H₂)' in name:
        color_mapping[surf_id] = shades['orange'][base_colors['orange']]
        base_colors['orange'] += 1
    elif 'vac' in name:
        color_mapping[surf_id] = shades['dodgerblue'][base_colors['dodgerblue']]
        base_colors['dodgerblue'] += 1
    elif 'clean' in name:
        color_mapping[surf_id] = shades['limegreen'][base_colors['limegreen']]
        base_colors['limegreen'] += 1
    elif '*OH+*OH' in name:
        color_mapping[surf_id] = shades['darkgoldenrod'][base_colors['darkgoldenrod']]
        base_colors['darkgoldenrod'] += 1
    elif '*OH+*O' in name:
        color_mapping[surf_id] = shades['mediumblue'][base_colors['mediumblue']]
        base_colors['mediumblue'] += 1
    elif '*O+*O' in name:
        color_mapping[surf_id] = shades['silver'][base_colors['silver']]
        base_colors['silver'] += 1
    elif '*OH' in name:
        color_mapping[surf_id] = shades['gold'][base_colors['gold']]
        base_colors['gold'] += 1
    elif '*O' in name:
        color_mapping[surf_id] = shades['hotpink'][base_colors['hotpink']]
        base_colors['hotpink'] += 1
    else:
        color_mapping[surf_id] = 'white'  # fallback color

# Apply color mapping to the colormap and ID mapping
colors = [color_mapping[sid] for sid in unique_ids]
cmap = mcolors.ListedColormap(colors)
bounds = np.arange(len(colors) + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Regenerate the mapped surface data
id_map = {val: idx for idx, val in enumerate(unique_ids)}
mapped_surfaces = np.vectorize(id_map.get)(lowest_surfaces)

# Create legend using the same surface ID mapping order
for idx, surf_id in enumerate(unique_ids):
    label = surfs[int(surf_id)][10]
    plt.plot([], [], color=colors[idx], linewidth=5, label=label)

# pcolormesh
pH, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH, U, mapped_surfaces, cmap=cmap, norm=norm)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., 
           fontsize='small', ncol=1, handlelength=3, edgecolor='black') 

plt.plot(pHrange, 1.23-pHrange*const, '--', color='black', lw=1)
ax.text(0.2, 1.00, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$', rotation=-12)
plt.plot(pHrange, 0-pHrange*const, '--', color='black', lw=1)
ax.text(0.2, -0.15 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$', rotation=-12)

plt.tight_layout()
if GCDFT:
    plt.savefig(f'pourbaixGC.png', dpi=300, bbox_inches='tight')
else:
    plt.savefig(f'pourbaix.png', dpi=300, bbox_inches='tight')
plt.show()