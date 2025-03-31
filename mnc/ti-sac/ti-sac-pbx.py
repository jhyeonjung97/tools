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
gas_nitrogen = -16.64503942/2 # from N₂, eV

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--gc', action='store_true', help='Enable GCDFT mode')
parser.add_argument('--bulk', action='store_true', help='Enable bulk Pourbaix mode')
parser.add_argument('--suffix', type=str, default='', help='Suffix for output filename')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--save-dir', action='store_true', help='Save to predefined directory')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
args = parser.parse_args()

GCDFT = args.gc
BULK_PB = args.bulk
target_pH = args.ph

if BULK_PB:
    ref_surf_name = 'vac+Ti(s)+N₂(g)/2'
else:
    ref_surf_name = '*'

if args.suffix:
    suffix = '_' + args.suffix
else:
    suffix = ''

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
pHrange = np.arange(0, 14+tick, tick)
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
        surface_term = ((surfs[k][8]*(U**2) + surfs[k][9]*U + surfs[k][10])
                          - (surfs[n_ref][8]*(U**2) + surfs[n_ref][9]*U + surfs[n_ref][10]))
    else:
        surface_term = surfs[k][0] - surfs[n_ref][0]
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
    [ -19.100/calmol, 0, 1, +0, 1, 0, 3, 0, 0, 0, 0, 'HNO₃(aq)'],
    [ -63.050/calmol, 0, 1, +0, 1, 1, 0, 0, 0, 0, 0, 'NHOH(aq)'],
    [ -19.000/calmol, 0, 1, +1, 1, 0, 0, 0, 0, 0, 0, 'NH⁺(aq)'],
    [  30.560/calmol, 0, 2, +0, 4, 0, 0, 0, 0, 0, 0, 'N₂H₄(aq)'],
    [  22.500/calmol, 0, 2, +2, 6, 0, 0, 0, 0, 0, 0, 'N₂H₆²⁺(aq)'],
    [  -5.600/calmol, 0, 1, +0, 2, 1, 0, 0, 0, 0, 0, 'NH₂OH(aq)'],
    [ -13.540/calmol, 0, 1, +0, 3, 1, 0, 0, 0, 0, 0, 'NH₂OH₂(aq)'],
    [  71.300/calmol, 0, 3, +0, 1, 0, 0, 0, 0, 0, 0, 'HN₃(aq)'],
    [  77.700/calmol, 0, 3, -1, 0, 0, 0, 0, 0, 0, 0, 'N₃⁻(aq)'],
    [   2.994/calmol, 0, 2, +0, 0, 0, 0, 0, 0, 0, 0, 'N₂(aq)'],
    [   8.600/calmol, 0, 2, +1, 2, 0, 2, 0, 0, 0, 0, 'H₂N₂O₂⁺(aq)'],
    [  33.200/calmol, 0, 2, -2, 0, 0, 2, 0, 0, 0, 0, 'N₂O₂²⁻(aq)'],
    [  18.200/calmol, 0, 1, -1, 2, 0, 2, 0, 0, 0, 0, 'NH₂O₂⁻(aq)'],
    [ -12.820/calmol, 0, 1, +0, 1, 0, 2, 0, 0, 0, 0, 'HNO₂(aq)'],
    [  -8.250/calmol, 0, 1, -1, 0, 0, 2, 0, 0, 0, 0, 'NO₂⁻(aq)'],
    [ -26.430/calmol, 0, 1, +0, 1, 0, 3, 0, 0, 0, 0, 'HNO₃(aq)'],
    [ -26.430/calmol, 0, 1, -1, 0, 0, 3, 0, 0, 0, 0, 'NO₃⁻(aq)'],
]

solids = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0,               1, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'Ti(s)'],
    [-116.920/calmol, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, 'TiO(s)'],
    [-342.310/calmol, 2, 0, +0, 0, 0, 3, 0, 0, 0, 0, 'Ti₂O₃(s)'],
    [-250.905/calmol, 1, 0, +0, 0, 3, 0, 0, 0, 0, 0, 'Ti(OH)₃(s)'],
    [-553.120/calmol, 3, 0, +0, 0, 0, 5, 0, 0, 0, 0, 'Ti₃O₅(s)'],
    [-212.330/calmol, 1, 0, +0, 0, 0, 2, 0, 0, 0, 0, 'TiO₂(s)'],
    [-252.990/calmol, 1, 0, +0, 1, 1, 2, 0, 0, 0, 0, 'TiO₂·H₂O(s)'],
    [  32.000/calmol, 0, 2, +0, 0, 0, 5, 0, 0, 0, 0, 'N₂O₅(l)'],
]

gases = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0            , 0, 2, +0, 0, 0, 0, 0, 0, 0, 0, 'N₂(g)'],
    [-3.976/calmol, 0, 1, +0, 3, 0, 0, 0, 0, 0, 0, 'NH₃(g)'],
    [78.500/calmol, 0, 3, +0, 1, 0, 0, 0, 0, 0, 0, 'HN₃(g)'],
    [81.476/calmol, 0, 1, +0, 0, 0, 0, 0, 0, 0, 0, 'N(g)'],
    [24.760/calmol, 0, 2, +0, 0, 0, 1, 0, 0, 0, 0, 'N₂O(g)'],
    [20.719/calmol, 0, 1, +0, 0, 0, 1, 0, 0, 0, 0, 'NO(g)'],
    [12.390/calmol, 0, 1, +0, 0, 0, 2, 0, 0, 0, 0, 'NO₂(g)'],
    [23.491/calmol, 0, 2, +0, 0, 0, 4, 0, 0, 0, 0, 'N₂O₄(g)'],
]

nions, nsolids, ngases = len(ions), len(solids), len(gases)
for i in range(nions):
    ions[i][0] += water * (ions[i][5] + ions[i][6] + 2*ions[i][7]) + bulk_metal * ions[i][1] + gas_nitrogen * ions[i][2]
    if ions[i][1] > 1 and ions[i][2] > 1:
        raise ValueError(f"Ion {ions[i][11]} has both Ti and N")
    if ions[i][1] > 1:
        ions[i][11] = ions[i][11] + f'/{int(ions[i][1])}'
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]
    elif ions[i][2] > 1:
        ions[i][11] = ions[i][11] + f'/{int(ions[i][2])}'
        ions[i] = [x / ions[i][2] if isinstance(x, (int, float)) else x for x in ions[i]]
    ions[i][10] = ions[i][0]

for s in range(nsolids):
    solids[s][0] += water * (solids[s][5] + solids[s][6] + 2*solids[s][7]) + bulk_metal * solids[s][1] + gas_nitrogen * solids[s][2]
    if solids[s][1] > 1 and solids[s][2] > 1:
        raise ValueError(f"Solid {solids[s][11]} has both Ti and N")
    if solids[s][1] > 1:
        solids[s][11] = solids[s][11] + f'/{int(solids[s][1])}'
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]
    elif solids[s][2] > 1:
        solids[s][11] = solids[s][11] + f'/{int(solids[s][2])}'
        solids[s] = [x / solids[s][2] if isinstance(x, (int, float)) else x for x in solids[s]]
    solids[s][10] = solids[s][0]

for g in range(ngases):
    gases[g][0] += water * (gases[g][5] + gases[g][6] + 2*gases[g][7]) + bulk_metal * gases[g][1] + gas_nitrogen * gases[g][2]
    if gases[g][1] > 1 and gases[g][2] > 1:
        raise ValueError(f"Gas {gases[g][11]} has both Ti and N")
    if gases[g][1] > 1:
        gases[g][11] = gases[g][11] + f'/{int(gases[g][1])}'
        gases[g] = [x / gases[g][1] if isinstance(x, (int, float)) else x for x in gases[g]]
    elif gases[g][2] > 1:
        gases[g][11] = gases[g][11] + f'/{int(gases[g][2])}'
        gases[g] = [x / gases[g][2] if isinstance(x, (int, float)) else x for x in gases[g]]
    gases[g][10] = gases[g][0]

surfs = [
    # ['E', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-269.569746,   0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'vac'],
    [-279.340765,   0, 0, +0, 2, 0, 0, 0, 0, 0, 0, 'vac(H₂)'],
    [-277.23249688, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*'],
    [-305.01556231, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(mono)'],
    [-305.14584901, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(bi)'],
    [-300.39230369, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(anti)'],
    [-301.21590570, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(adj)'],
    [-297.04551538, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*HNO'],
    [-310.30495026, 1, 1, +0, 0, 1, 2, 0, 0, 0, 0, '*NO₂(mono)+*OH(adj)'],
    [-300.07883828, 1, 1, +0, 1, 1, 0, 0, 0, 0, 0, '*HNOH'],
    [-283.58401076, 1, 0, +0, 2, 0, 0, 0, 0, 0, 0, '*H+*H(anti)'],
    [-284.73355007, 1, 0, +0, 2, 0, 0, 0, 0, 0, 0, '*H+*H(adj)'],
    [-294.10445417, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*H+*NO(anti)'],
    [-281.52877857, 1, 0, +0, 1, 0, 0, 0, 0, 0, 0, '*H'],
    [-284.81503764, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*N'],
    [-290.82616987, 1, 1, +0, 1, 0, 0, 0, 0, 0, 0, '*NH'],
    [-295.02930066, 1, 1, +0, 2, 0, 0, 0, 0, 0, 0, '*NH₂'],
    [-306.36204655, 1, 1, +0, 2, 1, 0, 0, 0, 0, 0, '*NH₂+*OH(adj)'],
    [-292.62371142, 1, 1, +0, 0, 0, 1, 0, 0, 0, 0, '*NO'],
    [-299.05807974, 1, 1, +0, 0, 0, 2, 0, 0, 0, 0, '*NO₂'],
    [-307.08538159, 1, 1, +0, 1, 0, 3, 0, 0, 0, 0, '*NO₃(bi)+*H(anti)'],
    [-296.12169107, 1, 1, +0, 0, 1, 0, 0, 0, 0, 0, '*NOH'],
    [-303.48187042, 1, 1, +0, 0, 1, 1, 0, 0, 0, 0, '*NO+*OH(adj)'],
    [-299.63088262, 1, 1, +0, 0, 0, 2, 0, 0, 0, 0, '*NO₂(bi)'],
    [-291.16160779, 1, 1, +0, 0, 0, 1, 0, 0, 0, 0, '*ON'],
    [-301.40150027, 1, 1, +0, 2, 0, 1, 0, 0, 0, 0, '*ONH₂'],
    [-306.42159787, 1, 1, +0, 3, 0, 1, 0, 0, 0, 0, '*O+*NH₃(adj)'],
    [-303.06694618, 1, 1, +0, 1, 0, 2, 0, 0, 0, 0, '*HNO₂'],
    [-295.88614166, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*ONH'],
    [-286.71126084, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
    [-290.31253022, 1, 0, +0, 0, 0, 2, 0, 0, 0, 0, '*O+*O(anti)'],
]

nsurfs = len(surfs)
ref_surf = surfs[0][0]
ref_surf_gc = surfs[0][10]
for k in range(nsurfs):
    formation_energy_corr = (
        - surfs[k][4] * (gh - dgh) # H
        - surfs[k][5] * (goh - dgoh) # OH
        - surfs[k][6] * (go - dgo) # O
        - surfs[k][7] * (gooh - dgooh) # OOH
    )
    surfs[k][0] = surfs[k][0] - ref_surf + formation_energy_corr 
    surfs[k][10] = surfs[k][10] - ref_surf_gc + formation_energy_corr 

if BULK_PB:
    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][1] == 0:
            for i, ion in enumerate(ions):
                if ion[1] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + ion[j])
                    new_surf.append(surfs[k][11] + '+' + ion[11])
                    new_surfs.append(new_surf)
            for s, solid in enumerate(solids):
                if solid[1] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + solid[j])
                    new_surf.append(surfs[k][11] + '+' + solid[11])
                    new_surfs.append(new_surf)
            for g, gas in enumerate(gases):
                if gas[1] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + gas[j])
                    new_surf.append(surfs[k][11] + '+' + gas[11])
                    new_surfs.append(new_surf)
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[1] != 0]
    nsurfs = len(surfs)

    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][2] == 0:
            for i, ion in enumerate(ions):
                if ion[2] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + ion[j])
                    new_surf.append(surfs[k][11] + '+' + ion[11])
                    new_surfs.append(new_surf)
            for s, solid in enumerate(solids):
                if solid[2] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + solid[j])
                    new_surf.append(surfs[k][11] + '+' + solid[11])
                    new_surfs.append(new_surf)
            for g, gas in enumerate(gases):
                if gas[2] == 1:
                    new_surf = []
                    for j in range(11):
                        new_surf.append(surfs[k][j] + gas[j])
                    new_surf.append(surfs[k][11] + '+' + gas[11])
                    new_surfs.append(new_surf)
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[2] != 0]
    nsurfs = len(surfs)

for i in range(nsurfs):
    if surfs[i][11] == ref_surf_name:
        n_ref = i
    print(f"#{i+1}:\t{surfs[i][0]:.2f}\t{surfs[i][11]}")

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

print("\nlowest_surfaces:")
print(lowest_surfaces)

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
    name = surfs[int(sid)][11]
    print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}, name = {name}")

# Set Axes Limits and Labels
fig, ax = plt.subplots(figsize=(7, 6), dpi=100)
ax.axis([0, 14, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

# Check unique values in lowest_surfaces and create a list of unique surface IDs
unique_ids = np.unique(lowest_surfaces)

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

print(f"\nSurfaces at pH={target_pH}:")
for sid in unique_ids_pH:
    name = surfs[int(sid)][11]
    print(f"Surface {sid}: {name}")

# Count surfaces for each color group
color_counts = {
    'RdYlBu': 0,    # vac 포함된 모든 표면
    'gray': 0,      # others
}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][11]
    if 'vac' in name:  # vac(H₂)도 포함
        color_counts['RdYlBu'] += 1
    else:
        color_counts['gray'] += 1

# Define base color groups and their initial indices
base_colors = {
    'RdYlBu': 0,    # vac 포함된 모든 표면
    'gray': 0,      # others
}

# Generate custom colormaps and shades
cmaps = {}
shades = {}

# RdYlBu는 matplotlib의 내장 colormap 사용
rdylbu_cmap = plt.get_cmap('RdYlBu')
if color_counts['RdYlBu'] > 0:
    shades['RdYlBu'] = [rdylbu_cmap(i) for i in np.linspace(0, 1, color_counts['RdYlBu'])]

# 나머지 색상은 기존 방식대로
for base_color in ['gray']:
    cmap = LinearSegmentedColormap.from_list(f"custom_{base_color}", ["white", base_color])
    cmaps[base_color] = cmap
    if color_counts[base_color] > 0:
        shades[base_color] = [cmap(i) for i in np.linspace(0.1, 0.9, color_counts[base_color])]
        
# Map surface ID to corresponding color shade
color_mapping = {}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][11]
    if 'vac' in name:  # vac(H₂)도 포함
        color_mapping[surf_id] = shades['RdYlBu'][base_colors['RdYlBu']]
        base_colors['RdYlBu'] += 1
    else:
        color_mapping[surf_id] = shades['gray'][base_colors['gray']]
        base_colors['gray'] += 1
        
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
    label = surfs[int(surf_id)][11]
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
png_name = 'pourbaix'
if BULK_PB:
    png_name += '_bulk'
else:
    png_name += '_surf'
if GCDFT:
    png_name += '_gc'

# Add potential vs. energy plot at pH=0
fig2, ax2 = plt.subplots(figsize=(7, 6), dpi=100)
ax2.axis([Umin, Umax, None, None])  # Fix x-axis only, y-axis will be adjusted automatically
ax2.set_xlabel('E (V vs. SHE)', labelpad=0)
ax2.set_ylabel('ΔG (eV)', labelpad=-6)
ax2.tick_params(right=True, direction="in")

# Calculate energy at specific pH
all_energies = []
unique_ids_set = set(unique_ids_pH)  # unique_ids를 set으로 변환하여 검색 속도 향상

for k in range(2, len(surfs)):
    energies = np.zeros(len(Urange))
    for i, U in enumerate(Urange):
        energy = dg(k, target_pH, U, concentration=1e-6, n_ref=n_ref)
        energies[i] = energy
    all_energies.extend(energies)
    
    if k in unique_ids_set:
        ax2.plot(Urange, energies, label=surfs[k][11], lw=0.5)
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

# Save original Pourbaix diagram
plt.figure(1)
plt.savefig(f'{save_dir}{png_name}{suffix}.png', dpi=300, bbox_inches='tight')
print(f"Saved plots to {png_name}{suffix}.png and {png_name}_pH{target_pH}{suffix}.png")

if args.show:
    plt.show()