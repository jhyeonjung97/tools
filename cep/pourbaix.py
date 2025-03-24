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
bulk_metal = -5.041720865 # Fe, eV

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
        surface_term = ((surfs[k][7]*(U**2) + surfs[k][8]*U + surfs[k][9])
                          - (surfs[0][7]*(U**2) + surfs[0][8]*U + surfs[0][9]))
    else:
        surface_term = surfs[k][0] - surfs[0][0]
    U_coeff = 1*surfs[k][3] - 1*surfs[k][4] - 2*surfs[k][5] - 3*surfs[k][6] - surfs[k][2]
    pH_coeff = 1*surfs[k][3] - 1*surfs[k][4] - 2*surfs[k][5] - 3*surfs[k][6]
    dg = surface_term + U_coeff * U + pH_coeff * const * pH
    if '(aq)' in surfs[k][10]:
        dg += const * log10(concentration)
    return dg
    
ions = [
    # ['Ef', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-20.300/calmol,  1, +2, 0, 0, 0, 0, 0, 0, 0, 'Fe²⁺(aq)'],
    [-90.627/calmol,  1, -1, 1, 0, 2, 0, 0, 0, 0, 'HFeO₂⁻(aq)'],
    [-2.530/calmol,   1, +3, 0, 0, 0, 0, 0, 0, 0, 'Fe³⁺(aq)'],
    [-55.910/calmol,  1, +2, 0, 1, 0, 0, 0, 0, 0, 'FeOH²⁺(aq)'],
    [-106.200/calmol, 1, +1, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂⁺(aq)'],
]

solids = [
    # ['Ef', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0,               1, +0, 0, 0, 0, 0, 0, 0, 0, 'Fe(s)'],
    [-58.880/calmol,  1, +0, 0, 0, 1, 0, 0, 0, 0, 'FeO'],
    [-242.400/calmol, 3, +0, 0, 0, 4, 0, 0, 0, 0, 'Fe₃O₄'],
    # [-177.100/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-161.930/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-115.570/calmol, 1, +0, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂'],
    [-166.000/calmol, 1, +0, 0, 3, 0, 0, 0, 0, 0, 'Fe(OH)₃'],
]

nions, nsolids = len(ions), len(solids)
for i in range(nions):
    ions[i][0] += water * (ions[i][4] + ions[i][5] + 2*ions[i][6]) + bulk_metal * ions[i][1]
    if ions[i][1] > 1:
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]

for s in range(nsolids):
    solids[s][0] += water * (solids[s][4] + solids[s][5] + 2*solids[s][6]) + bulk_metal * solids[s][1]
    if solids[s][1] > 1:
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]

surfs = [
    # ['E', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-269.569746, 0, +0, 0, 0, 0, 0, -0.3655405668135691, 0.45552073788487135, -269.645443005364, 'vac'],
    [-279.340765, 0, +0, 2, 0, 0, 0, -0.4310220085501212, -0.18813647315805954, -279.27741924851745, 'vac(H₂)'],
    [-277.215281, 1, +0, 0, 0, 0, 0, -0.5088279832115933, -0.3089382644866646, -277.16659314038986, 'clean(HS)'],
    [-277.546212, 1, +0, 0, 0, 0, 0, -0.423357238170348, -0.18029803360057767, -277.46155128048673, 'clean(IS)'],
    [-275.904061, 1, +0, 0, 0, 0, 0, -0.423198992129675, -0.023209968264605602, -275.8141779964752, 'clean(LS)'],
    [-288.221008, 1, +0, 0, 1, 0, 0, -0.4467206439644741, 0.04059102541473414, -288.1146816835875, '*OH(HS)'],
    [-288.124300, 1, +0, 0, 1, 0, 0, -0.49213990645417277, 0.2598149164035292, -288.04840908267477, '*OH(IS)'],
    [-287.563089, 1, +0, 0, 1, 0, 0, -0.4807076337804513, 0.18746381957939667, -287.4849251837713, '*OH(LS)'],
    [-283.228986, 1, +0, 0, 0, 1, 0, -0.40383700613414103, 0.2916884347090441, -283.1526637494537, '*O(HS)'],
    [-282.942139, 1, +0, 0, 0, 1, 0, -0.45657630964544965, 0.3430163672241705, -282.89037791515415, '*O(IS)'],
    [-282.513697, 1, +0, 0, 0, 1, 0, -0.4407271111295996, 0.2676105696973955, -282.43318248088525, '*O(LS)'],
    [-297.295568, 1, +0, 0, 2, 0, 0, -0.5595101320515073, 0.4931725265118078, -297.292789623116, '*OH+*OH(HS)'],
    [-296.997610, 1, +0, 0, 2, 0, 0, -0.5277748456986603, 0.39666069296217893, -297.15342439890617, '*OH+*OH(IS)'],
    [-296.898685, 1, +0, 0, 2, 0, 0, -0.59572540226707, 0.5697551318500524, -297.0414929981233, '*OH+*OH(LS)'],
    [-292.18643175, 1, +0, 0, 1, 1, 0, -0.474403700041583, 0.290084262601454, -292.08048960254666, '*OH+*O(HS)'],
    [-292.66930911, 1, +0, 0, 1, 1, 0, -0.46048018108346694, 0.6326720159755281, -292.7217143067089, '*OH+*O(LS)'],
    [-286.22738149, 1, +0, 0, 0, 2, 0, -0.5129915832564392, 0.7798114189274047, -286.3805242333587, '*O+*O(HS)'],
    [-286.76062327, 1, +0, 0, 0, 2, 0, -0.4870021793211474, 0.7400411428427354, -286.902324627179, '*O+*O(LS)'],
]

nsurfs = len(surfs)
surface_reference = surfs[0][0]
for k in range(nsurfs):
    formation_energy_corr = (
        - surface_reference
        - surfs[k][3] * (gh - dgh) # H
        - surfs[k][4] * (goh - dgoh) # OH
        - surfs[k][5] * (go - dgo) # O
        - surfs[k][6] * (gooh - dgooh) # OOH
    )
    surfs[k][0] += formation_energy_corr # E
    surfs[k][9] += formation_energy_corr # C
    
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
fig, ax = plt.subplots(figsize=(8, 5), dpi=100)
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

plt.plot(pHrange, 1.23-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 1.00, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$',
        color='blue', rotation=-12, fontsize=12)
plt.plot(pHrange, 0-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, -0.15 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$',
        color='blue', rotation=-12, fontsize=12)

plt.tight_layout()
if GCDFT:
    plt.savefig(f'pourbaixGC.png', dpi=300, bbox_inches='tight')
else:
    plt.savefig(f'pourbaix.png', dpi=300, bbox_inches='tight')
plt.show()