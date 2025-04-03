### Hailey ###
### Date: 2025-03-13 ###

import os
import numpy as np
from math import log10
from matplotlib import colormaps
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import argparse
import socket

bulk_metal = -5.041720865 # Fe, eV

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--gc', action='store_true', help='Enable GCDFT mode')
parser.add_argument('--bulk', action='store_true', help='Enable bulk Pourbaix mode')
parser.add_argument('--suffix', type=str, default='', help='Suffix for output filename')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--save-dir', action='store_true', help='Save to predefined directory')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--x', type=float, default=7, help='Figure width in inches (default: 7)')
parser.add_argument('--y', type=float, default=6, help='Figure height in inches (default: 6)')
args = parser.parse_args()

GCDFT = args.gc
BULK_PB = args.bulk
target_pH = args.ph

if args.suffix:
    suffix = '_' + args.suffix
else:
    suffix = ''

png_name = 'pourbaix'
if BULK_PB:
    png_name += '_bulk'
else:
    png_name += '_surf'
if GCDFT:
    png_name += '_gc'
    
if not args.save_dir:
    save_dir = ''
else:
    hostname = socket.gethostname()
    user_name = os.getlogin()
    if hostname == 'PC102616':
        save_dir = '/Users/jiuy97/Desktop/9_pourbaixGC/figures/'
    elif user_name == 'jiuy97':
        save_dir = '/pscratch/sd/j/jiuy97/9_pourbaixGC/figures/'
    elif user_name == 'hailey':
        save_dir = '/Users/hailey/Desktop/9_pourbaixGC/figures/'
    else:
        raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")
        
# constants
kb = 8.617e-5 # eV/K
T = 298.15 # K
const = kb * T * np.log(10) # 0.0592 eV
water = 2.4583 # the standard Gibbs free energy of formation of water

# units
kjmol = 96.485
calmol = 23.061

# ticks
tick = args.tick
pHmin, pHmax = 0, 14
pHrange = np.arange(pHmin, pHmax + tick, tick)
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

def dg(k, pH, U, concentration, n_ref):
    if GCDFT:
        surface_term = ((surfs[k][7]*(U**2) + surfs[k][8]*U + surfs[k][9])
                          - (surfs[n_ref][7]*(U**2) + surfs[n_ref][8]*U + surfs[n_ref][9]))
    else:
        surface_term = surfs[k][0] - surfs[n_ref][0]
    U_term = 1*surfs[k][3] -1*surfs[k][4] -2*surfs[k][5] -3*surfs[k][6] -surfs[k][2]
    pH_term = 1*surfs[k][3] -1*surfs[k][4] -2*surfs[k][5] -3*surfs[k][6]
    dg = surface_term + U_term * U + pH_term * const * pH
    if '(aq)' in surfs[k][10]:
        dg += const * log10(concentration)
    return dg
    
ions = [
    # ['Ef', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [ -20.300/calmol, 1, +2, 0, 0, 0, 0, 0, 0, 0, 'Fe²⁺'],
    [ -90.627/calmol, 1, -1, 1, 0, 2, 0, 0, 0, 0, 'HFeO₂⁻'],
    [  -2.530/calmol, 1, +3, 0, 0, 0, 0, 0, 0, 0, 'Fe³⁺'],
    [ -55.910/calmol, 1, +2, 0, 1, 0, 0, 0, 0, 0, 'FeOH²⁺'],
    [-106.200/calmol, 1, +1, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂⁺'],
    [ -40.000/calmol, 1, -1, 0, 0, 4, 0, 0, 0, 0, 'FeO₄⁻'],
]

solids = [
    # ['Ef', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0,               1, +0, 0, 0, 0, 0, 0, 0, 0, 'Fe'],
    [-58.880/calmol,  1, +0, 0, 0, 1, 0, 0, 0, 0, 'FeO'],
    [-242.400/calmol, 3, +0, 0, 0, 4, 0, 0, 0, 0, 'Fe₃O₄'],
    [-177.100/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-161.930/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-115.570/calmol, 1, +0, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂'],
    [-166.000/calmol, 1, +0, 0, 3, 0, 0, 0, 0, 0, 'Fe(OH)₃'],
]

nions, nsolids = len(ions), len(solids)
for i in range(nions):
    ions[i][0] += water * (ions[i][4] + ions[i][5] + 2*ions[i][6]) + bulk_metal * ions[i][1]
    if ions[i][1] > 1:
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]
    ions[i][9] = ions[i][0]
    ions[i][10] = '+' + ions[i][10] + '(aq)'

for s in range(nsolids):
    solids[s][0] += water * (solids[s][4] + solids[s][5] + 2*solids[s][6]) + bulk_metal * solids[s][1]
    if solids[s][1] > 1:
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]
    solids[s][9] = solids[s][0]
    solids[s][10] = '+' + solids[s][10] + '(s)'

# surfs = [
#     # ['E', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
#     [-271.953170, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'vac'],
#     [-281.836091, 0, +0, 2, 0, 0, 0, 0, 0, 0, 'vac(H₂)'],
#     [-280.176972, 1, +0, 0, 0, 0, 0, 0, 0, 0, 'clean'],
#     [-282.637740, 1, +0, 1, 0, 0, 0, 0, 0, 0, '*H'],
#     [-290.947634, 1, +0, 0, 1, 0, 0, 0, 0, 0, '*OH'],
#     [-285.944778, 1, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
#     [-300.626405, 1, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH'],
#     [-295.526586, 1, +0, 0, 1, 1, 0, 0, 0, 0, '*OH+*O'],
#     [-289.676342, 1, +0, 0, 0, 2, 0, 0, 0, 0, '*O+*O'],
# ]

surfs = [
    # ['E', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-269.569746, 0, +0, 0, 0, 0, 0, -0.365540566813569,  0.4555207378848713, -269.645443005364, 'vac'],
    [-279.340765, 0, +0, 2, 0, 0, 0, -0.431022008550121, -0.1881364731580595, -279.277419248517, 'vac(H₂)'],
    [-277.215281, 1, +0, 0, 0, 0, 0, -0.508827983211593, -0.3089382644866646, -277.166593140389, 'clean(HS)'],
    [-277.546212, 1, +0, 0, 0, 0, 0, -0.423357238170348, -0.1802980336005776, -277.461551280486, 'clean(IS)'],
    [-275.904061, 1, +0, 0, 0, 0, 0, -0.423198992129675, -0.0232099682646056, -275.814177996475, 'clean(LS)'],
    [-288.221008, 1, +0, 0, 1, 0, 0, -0.446720643964474,  0.0405910254147341, -288.114681683587, '*OH(HS)'],
    [-288.124300, 1, +0, 0, 1, 0, 0, -0.492139906454172,  0.2598149164035292, -288.048409082674, '*OH(IS)'],
    [-287.563089, 1, +0, 0, 1, 0, 0, -0.480707633780451,  0.1874638195793966, -287.484925183771, '*OH(LS)'],
    [-283.228986, 1, +0, 0, 0, 1, 0, -0.403837006134141,  0.2916884347090441, -283.152663749453, '*O(HS)'],
    [-282.942139, 1, +0, 0, 0, 1, 0, -0.456576309645449,  0.3430163672241705, -282.890377915154, '*O(IS)'],
    [-282.513697, 1, +0, 0, 0, 1, 0, -0.440727111129599,  0.2676105696973955, -282.433182480885, '*O(LS)'],
    [-297.295568, 1, +0, 0, 2, 0, 0, -0.559510132051507,  0.4931725265118078, -297.292789623116, '*OH+*OH(HS)'],
    [-296.997610, 1, +0, 0, 2, 0, 0, -0.527774845698660,  0.3966606929621789, -297.153424398906, '*OH+*OH(IS)'],
    [-296.898685, 1, +0, 0, 2, 0, 0, -0.595725402267070,  0.5697551318500524, -297.041492998123, '*OH+*OH(LS)'],
    [-292.186431, 1, +0, 0, 1, 1, 0, -0.474403700041583,  0.2900842626014540, -292.080489602546, '*OH+*O(HS)'],
    [-292.669309, 1, +0, 0, 1, 1, 0, -0.460480181083466,  0.6326720159755281, -292.721714306708, '*OH+*O(LS)'],
    [-286.227381, 1, +0, 0, 0, 2, 0, -0.512991583256439,  0.7798114189274047, -286.380524233358, '*O+*O(HS)'],
    [-286.760623, 1, +0, 0, 0, 2, 0, -0.487002179321147,  0.7400411428427354, -286.902324627179, '*O+*O(LS)'],
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

if BULK_PB:
    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][1] == 0:
            for i in range(nions):
                new_surf = []
                for j in range(11):
                    new_surf.append(surfs[k][j] + ions[i][j])
                new_surfs.append(new_surf)
            for s in range(nsolids):
                new_surf = []
                for j in range(11):
                    new_surf.append(surfs[k][j] + solids[s][j])
                new_surfs.append(new_surf)
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[1] != 0]
    nsurfs = len(surfs)
else:
    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][1] == 0:
            new_surf = []
            for j in range(11):
                new_surf.append(surfs[k][j] + solids[0][j])
            new_surfs.append(new_surf)    
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[1] != 0]
    nsurfs = len(surfs)

if GCDFT:
    print(f"No.\tEnergy\t#Fe\t#e\t#H\t#OH\t#O\tA\tB\tC\tSurface")
else:
    print(f"No.\tEnergy\t#Fe\t#e\t#H\t#OH\t#O\tSurface")
for i in range(nsurfs):
    if surfs[i][10] == 'vac'+solids[0][10]:
        n_ref = i
    if GCDFT:
        print(f"#{i+1}:\t{surfs[i][0]:.2f}\t{surfs[i][1]:.2f}\t{surfs[i][2]:.2f}\t{surfs[i][3]:.2f}\t{surfs[i][4]:.2f}\t{surfs[i][5]:.2f}\t{surfs[i][7]:.2f}\t{surfs[i][8]:.2f}\t{surfs[i][9]:.2f}\t{surfs[i][10]}")
    else:
        print(f"#{i+1}:\t{surfs[i][0]:.2f}\t{surfs[i][1]:.2f}\t{surfs[i][2]:.2f}\t{surfs[i][3]:.2f}\t{surfs[i][4]:.2f}\t{surfs[i][5]:.2f}\t{surfs[i][10]}")

print(f"total surfaces: {nsurfs}")
print(f"reference surface: {n_ref+1}")
lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

pHindex = 0
for pH in pHrange:
    Uindex = 0
    for U in Urange:
        values = []
        for k in range(nsurfs):
            value = dg(k, pH, U, concentration=1e-6, n_ref=n_ref)
            values.append(value)
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Uindex][pHindex] = sorted_values[0]
        Uindex+=1
    pHindex+=1

min_coords = {}
n_rows, n_cols = lowest_surfaces.shape

for j in range(n_cols):
    for i in range(n_rows):
        sid = int(lowest_surfaces[i, j])
        x = pHrange[j]
        y = Urange[i]
        if sid not in min_coords:
            min_coords[sid] = (x, y)
        else:
            current_x, current_y = min_coords[sid]
            if x < current_x or (x == current_x and y < current_y):
                min_coords[sid] = (x, y)

for sid in sorted(min_coords):
    x, y = min_coords[sid]
    name = surfs[int(sid)][10]
    print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}, name = {name}")

# Set Axes Limits and Labels
fig, ax = plt.subplots(figsize=(args.x, args.y), dpi=100)
ax.axis([pHmin, pHmax, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

# Check unique values in lowest_surfaces and create a list of unique surface IDs
unique_ids = np.unique(lowest_surfaces)

# Count surfaces for each color group
color_counts = {
    'YlOrBr': 0,    # surfaces containing 'vac'
    'Greys': 0,     # others
}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][10]
    if 'vac' in name:  # including 'vac(H₂)'
        color_counts['YlOrBr'] += 1
    else:
        color_counts['Greys'] += 1

# Define base color groups and their initial indices
base_colors = {
    'YlOrBr': 0,    # surfaces containing 'vac'
    'Greys': 0,     # others
}

# Generate custom colormaps and shades
cmaps = {}
shades = {}

# Use matplotlib's built-in colormaps for YlOrBr and Greys
ylorbr_cmap = plt.get_cmap('YlOrBr')
if color_counts['YlOrBr'] > 0:
    shades['YlOrBr'] = [ylorbr_cmap(i) for i in np.linspace(0.1, 0.5, color_counts['YlOrBr'])]

greys_cmap = plt.get_cmap('Greys')
if color_counts['Greys'] > 0:
    shades['Greys'] = [greys_cmap(i) for i in np.linspace(0.1, 0.5, color_counts['Greys'])]
        
# Map surface ID to corresponding color shade
color_mapping = {}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][10]
    if 'vac' in name:  # including 'vac(H₂)'
        color_mapping[surf_id] = shades['YlOrBr'][base_colors['YlOrBr']]
        base_colors['YlOrBr'] += 1
    else:
        color_mapping[surf_id] = shades['Greys'][base_colors['Greys']]
        base_colors['Greys'] += 1
        
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
pH_grid, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH_grid, U, mapped_surfaces, cmap=cmap, norm=norm)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., 
           fontsize='small', ncol=1, handlelength=3, edgecolor='black') 

plt.plot(pHrange, 1.23-pHrange*const, '--', lw=1, color='mediumblue')
ax.text(0.2, 1.23+0.12, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$', rotation=-13, color='mediumblue', ha='left', va='top')
plt.plot(pHrange, 0-pHrange*const, '--', lw=1, color='mediumblue')
ax.text(0.2, 0+0.12, r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$', rotation=-13, color='mediumblue', ha='left', va='top')

plt.tight_layout()
plt.savefig(f'{save_dir}{png_name}{suffix}.png', dpi=300, bbox_inches='tight')

# Add potential vs. energy plot at pH=0
fig2, ax2 = plt.subplots(figsize=(args.x, args.y), dpi=100)
ax2.axis([Umin, Umax, None, None])  # Fix x-axis only, y-axis will be adjusted automatically
ax2.set_xlabel('E (V vs. SHE)', labelpad=0)
ax2.set_ylabel('ΔG (eV)', labelpad=-6)
ax2.tick_params(right=True, direction="in")

# Calculate energy at specific pH
all_energies = []

# Calculate lowest surfaces at target_pH directly using dg function
lowest_surfaces_pH = np.zeros(len(Urange))
for Uindex, U in enumerate(Urange):
    values = []
    for k in range(nsurfs):
        value = dg(k, pH=target_pH, U=U, concentration=1e-6, n_ref=n_ref)
        values.append(value)
    sorted_values = sorted(range(len(values)), key=lambda k: values[k])
    lowest_surfaces_pH[Uindex] = sorted_values[0]

unique_ids_pH = np.unique(lowest_surfaces_pH.astype(int))
unique_ids_set = set(unique_ids_pH)  # unique_ids를 set으로 변환하여 검색 속도 향상

for k in range(nsurfs):
    energies = np.zeros(len(Urange))
    for i, U in enumerate(Urange):
        energy = dg(k, target_pH, U, concentration=1e-6, n_ref=n_ref)
        energies[i] = energy
    all_energies.extend(energies)
    
    if k in unique_ids_set:
        ax2.plot(Urange, energies, label=surfs[k][10], lw=0.5)
    else:
        ax2.plot(Urange, energies, lw=0.5)

# Adjust y-axis range
y_min = min(all_energies)
y_max = max(all_energies)
y_margin = (y_max - y_min) * 0.1
ax2.set_ylim(y_min - y_margin, y_max + y_margin)

ax2.set_xlim(Umin, Umax)
ax2.grid(True, linestyle='--', alpha=0.3)

ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., 
           fontsize='small', ncol=1, handlelength=3, edgecolor='black')

plt.tight_layout()
plt.savefig(f'{save_dir}{png_name}_pH{target_pH}{suffix}.png', dpi=300, bbox_inches='tight')
print(f"Saved plots to {png_name}{suffix}.png and {png_name}_pH{target_pH}{suffix}.png")

if args.show:
    plt.show()