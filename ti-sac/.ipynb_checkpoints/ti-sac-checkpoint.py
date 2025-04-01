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
parser.add_argument('--suffix', type=str, default='TiSAC0', help='Suffix for output filename')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--save-dir', action='store_true', help='Save to predefined directory')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
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
pHrange = np.arange(0, 14+tick, tick)
Umin, Umax = -2.0, 2.0
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
    dg = surface_term + U_coeff * U + (pH_coeff) * const * pH
    if '(aq)' in surfs[k][11] or '(g)' in surfs[k][11]:
        dg += const * log10(concentration)
    return dg
    
ions = [
    # ['Ef', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [ -75.100/calmol, 1, 0, +2, 0, 0, 0, 0, 0, 0, 0, 'Ti²⁺(aq)'],
    [ -83.600/calmol, 1, 0, +3, 0, 0, 0, 0, 0, 0, 0, 'Ti³⁺(aq)'],
    [-138.000/calmol, 1, 0, +2, 0, 0, 1, 0, 0, 0, 0, 'TiO²⁺(aq)'],
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
        ions[i][11] = f'{(1/ions[i][1]):.2f}' + ions[i][11]
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]
    elif ions[i][2] > 1:
        ions[i][11] = f'{(1/ions[i][2]):.2f}' + ions[i][11]
        ions[i] = [x / ions[i][2] if isinstance(x, (int, float)) else x for x in ions[i]]
    ions[i][10] = ions[i][0]

for s in range(nsolids):
    solids[s][0] += water * (solids[s][5] + solids[s][6] + 2*solids[s][7]) + bulk_metal * solids[s][1] + gas_nitrogen * solids[s][2]
    if solids[s][1] > 1 and solids[s][2] > 1:
        raise ValueError(f"Solid {solids[s][11]} has both Ti and N")
    if solids[s][1] > 1:
        solids[s][11] = f'{(1/solids[s][1]):.2f}' + solids[s][11]
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]
    elif solids[s][2] > 1:
        solids[s][11] = f'{(1/solids[s][2]):.2f}' + solids[s][11]
        solids[s] = [x / solids[s][2] if isinstance(x, (int, float)) else x for x in solids[s]]
    solids[s][10] = solids[s][0]

for g in range(ngases):
    gases[g][0] += water * (gases[g][5] + gases[g][6] + 2*gases[g][7]) + bulk_metal * gases[g][1] + gas_nitrogen * gases[g][2]
    if gases[g][1] > 1 and gases[g][2] > 1:
        raise ValueError(f"Gas {gases[g][11]} has both Ti and N")
    if gases[g][1] > 1:
        gases[g][11] = f'{(1/gases[g][1]):.2f}' + gases[g][11]
        gases[g] = [x / gases[g][1] if isinstance(x, (int, float)) else x for x in gases[g]]
    elif gases[g][2] > 1:
        gases[g][11] = f'{(1/gases[g][2]):.2f}' + gases[g][11]
        gases[g] = [x / gases[g][2] if isinstance(x, (int, float)) else x for x in gases[g]]
    gases[g][10] = gases[g][0]

surfs = [
    # ['E', '#M(Ti)', '#N', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-269.569746,   0, 0, +0, 0, 0, 0, 0, 0, 0, 0, 'vac'],
    [-279.340765,   0, 0, +0, 2, 0, 0, 0, 0, 0, 0, 'vac(H₂)'],
    # [-290.30625044, 1, 0, +0, 0, 1, 0, 0, 0, 0, 0, '*OH'],
    # [-286.71959099, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*O'],
    # [-301.53899925, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(adj)'],
    # [-300.39500297, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(anti)'],
    # [-295.43972963, 1, 0, +0, 0, 1, 1, 0, 0, 0, 0, '*O+*OH(anti)'],
    # [-295.43511415, 1, 0, +0, 0, 1, 1, 0, 0, 0, 0, '*OH+*O(anti)'],
    # [-289.47051131, 1, 0, +0, 0, 0, 2, 0, 0, 0, 0, '*O+*O(anti)'],
    # [-291.63620335, 1, 0, +0, 0, 0, 2, 0, 0, 0, 0, '*O₂'],
    # [-300.01076906, 1, 0, +0, 0, 0, 1, 1, 0, 0, 0, '*O+*OOH(anti)'],
    # [-304.62550385, 1, 0, +0, 0, 1, 0, 1, 0, 0, 0, '*OH+*OOH(anti)'],
    # [-294.51658115, 0, 0, +0, 0, 0, 0, 1, 0, 0, 0, '*OOH'],
    # [-309.56150674, 0, 0, +0, 0, 0, 0, 2, 0, 0, 0, '*OOH+*OOH(adj)'],
    [-304.97333512, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(mono)'],
    [-305.14580502, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(bi)'],
    [-277.22964386, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*'],
    [-277.23249688, 1, 0, +0, 0, 0, 0, 0, 0, 0, 0, '*'],
    [-292.89799477, 1, 0, +0, 0, 1, 1, 0, 0, 0, 0, '*OH₂'],
    [-281.52525144, 1, 0, +0, 1, 0, 0, 0, 0, 0, 0, '*H'],
    [-283.58400933, 1, 0, +0, 2, 0, 0, 0, 0, 0, 0, '*H+*H(anti)'],
    [-284.73355007, 1, 0, +0, 2, 0, 0, 0, 0, 0, 0, '*H+*H(adj)'],
    [-308.95411743, 1, 1, +0, 1, 0, 3, 0, 0, 0, 0, '*HNO₃(bi)'],
    [-294.10389720, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*NO+*H(anti)'],
    [-297.04728106, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*HNO'],
    [-284.81503764, 1, 1, +0, 1, 1, 0, 0, 0, 0, 0, '*HNOH'],
    [-300.07498499, 1, 1, +0, 0, 1, 1, 0, 0, 0, 0, '*HONO'],
    [-284.81503764, 1, 1, +0, 0, 0, 0, 0, 0, 0, 0, '*N'],
    [-295.02930066, 1, 1, +0, 2, 0, 0, 0, 0, 0, 0, '*NH₂'],
    [-306.36204655, 1, 1, +0, 2, 1, 0, 0, 0, 0, 0, '*NH₂+*OH(adj)'],
    [-298.07547708, 1, 1, +0, 3, 0, 0, 0, 0, 0, 0, '*NH₃'],
    [-310.32808382, 1, 1, +0, 3, 1, 0, 0, 0, 0, 0, '*NH₃+*OH(adj)'],
    [-299.05807974, 1, 1, +0, 0, 0, 2, 0, 0, 0, 0, '*NO₂(mono)'],
    [-305.80131325, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(bi)'],
    [-307.08538159, 1, 1, +0, 1, 0, 3, 0, 0, 0, 0, '*NO₃(bi)+*H(anti)'],
    [-292.63294558, 1, 1, +0, 0, 0, 1, 0, 0, 0, 0, '*NO'],
    [-296.12308988, 1, 1, +0, 0, 1, 0, 0, 0, 0, 0, '*NOH'],
    [-305.01556231, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(mono)'],
    [-286.71347573, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
    [-290.83585688, 1, 0, +0, 0, 1, 0, 0, 0, 0, 0, '*OH'],
    [-303.48187042, 1, 1, +0, 0, 1, 1, 0, 0, 0, 0, '*NO+*OH(adj)'],
    [-300.39230369, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(anti)'],
    [-301.21590570, 1, 0, +0, 0, 2, 0, 0, 0, 0, 0, '*OH+*OH(adj)'],
    [-281.52877857, 1, 0, +0, 1, 0, 0, 0, 0, 0, 0, '*H'],
    [-310.30495021, 1, 1, +0, 0, 1, 2, 0, 0, 0, 0, '*NO₂(mono)+*OH(adj)'],
    [-290.82616987, 1, 1, +0, 1, 0, 0, 0, 0, 0, 0, '*NH'],
    [-286.71345207, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
    [-290.31287362, 1, 0, +0, 0, 1, 0, 0, 0, 0, 0, '*OH'],
    [-305.24167142, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(mono)'],
    [-308.08067363, 1, 1, +0, 1, 0, 3, 0, 0, 0, 0, '*NO₃H(mono)'],
    [-299.63088262, 1, 1, +0, 0, 0, 2, 0, 0, 0, 0, '*NO₂(bi)'],
    [-305.24823013, 1, 1, +0, 0, 0, 3, 0, 0, 0, 0, '*NO₃(bi)'],
    [-308.47508401, 1, 1, +0, 1, 0, 3, 0, 0, 0, 0, '*NO₃H(bi)'],
    [-306.64575220, 1, 1, +0, 1, 1, 1, 0, 0, 0, 0, '*NHO₂H(bi)'],
    [-303.06694617, 1, 1, +0, 1, 0, 2, 0, 0, 0, 0, '*NHO₂(bi)'],
    [-302.62353824, 1, 1, +0, 0, 1, 1, 0, 0, 0, 0, '*NO₂H(bi)'],
    [-291.16160779, 1, 1, +0, 0, 0, 1, 0, 0, 0, 0, '*ON'],
    [-301.40150027, 1, 1, +0, 2, 0, 1, 0, 0, 0, 0, '*ONH₂'],
    # [-306.42159787, 1, 1, +0, 2, 0, 1, 0, 0, 0, 0, '*NH₃+*O(adj)'],
    [-295.88587495, 1, 1, +0, 1, 0, 1, 0, 0, 0, 0, '*ONH'],
    [-286.71126084, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
    [-286.47933739, 1, 0, +0, 0, 0, 1, 0, 0, 0, 0, '*O'],
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
else:
    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][1] == 0:
            new_surf = []
            for j in range(11):
                new_surf.append(surfs[k][j] + solids[0][j])
            new_surf.append(surfs[k][11] + '+' + solids[0][11])
            new_surfs.append(new_surf)
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[1] != 0]
    nsurfs = len(surfs)

    new_surfs = []
    for k in range(nsurfs):
        if surfs[k][2] == 0:
            new_surf = []
            for j in range(11):
                new_surf.append(surfs[k][j] + gases[0][j])
            new_surf.append(surfs[k][11] + '+' + gases[0][11])
            new_surfs.append(new_surf)
    surfs.extend(new_surfs)
    surfs = [surf for surf in surfs if surf[2] != 0]
    nsurfs = len(surfs)

print(f"No.\tEnergy\t#Fe\t#N\t#e\t#H\t#OH\t#O\tSurface")
for i in range(nsurfs):
    if surfs[i][11] == 'vac'+'+'+solids[0][11]+'+'+gases[0][11]:
        n_ref = i
    print(f"#{i+1}:\t{surfs[i][0]:.2f}\t{surfs[i][1]:.2f}\t{surfs[i][2]:.2f}\t{surfs[i][3]:.2f}\t{surfs[i][4]:.2f}\t{surfs[i][5]:.2f}\t{surfs[i][6]:.2f}\t{surfs[i][11]}")

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

# print("\nlowest_surfaces:")
# print(lowest_surfaces)

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
    'Blues': 0,     # surfaces containing 'vac'
    'Greys': 0,     # others
}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][11]
    if 'vac' in name:  # including 'vac(H₂)'
        color_counts['Blues'] += 1
    else:
        color_counts['Greys'] += 1

# Define base color groups and their initial indices
base_colors = {
    'Blues': 0,     # surfaces containing 'vac'
    'Greys': 0,     # others
}

# Generate custom colormaps and shades
cmaps = {}
shades = {}

# Use matplotlib's built-in colormaps for Blues and Greys
blues_cmap = plt.get_cmap('Blues')
if color_counts['Blues'] > 0:
    shades['Blues'] = [blues_cmap(i) for i in np.linspace(0.1, 0.5, color_counts['Blues'])]

greys_cmap = plt.get_cmap('Greys')
if color_counts['Greys'] > 0:
    shades['Greys'] = [greys_cmap(i) for i in np.linspace(0.1, 0.5, color_counts['Greys'])]
        
# Map surface ID to corresponding color shade
color_mapping = {}

for surf_id in unique_ids:
    name = surfs[int(surf_id)][11]
    if 'vac' in name:  # including 'vac(H₂)'
        color_mapping[surf_id] = shades['Blues'][base_colors['Blues']]
        base_colors['Blues'] += 1
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
    label = surfs[int(surf_id)][11]
    plt.plot([], [], color=colors[idx], linewidth=5, label=label)

# pcolormesh
pH_grid, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH_grid, U, mapped_surfaces, cmap=cmap, norm=norm)

plt.legend(bbox_to_anchor=(0.02, 1.02), loc='lower left', borderaxespad=0., 
           fontsize='small', ncol=2, handlelength=3, edgecolor='black') 

plt.plot(pHrange, 1.23-pHrange*const, '--', lw=1, color='mediumblue')
ax.text(0.2, 1.23+0.14, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$', rotation=-9, color='mediumblue', ha='left', va='top')
plt.plot(pHrange, 0-pHrange*const, '--', lw=1, color='mediumblue')
ax.text(0.2, 0+0.14, r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$', rotation=-9, color='mediumblue', ha='left', va='top')

plt.tight_layout()
plt.savefig(f'{save_dir}{png_name}{suffix}.png', dpi=300, bbox_inches='tight')

# Add potential vs. energy plot at pH=0
fig2, ax2 = plt.subplots(figsize=(7, 6), dpi=100)
ax2.axis([Umin, Umax, None, None])
ax2.set_xlabel('E (V vs. SHE)', labelpad=0)
ax2.set_ylabel('ΔG (eV)', labelpad=-6)
ax2.tick_params(right=True, direction="in")

# Calculate energy at specific pH
all_energies = []
unique_ids_set = set(unique_ids_pH)

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

ax2.legend(bbox_to_anchor=(0.02, 1.02), loc='lower left', borderaxespad=0., 
          fontsize='small', ncol=2, handlelength=3, edgecolor='black')

plt.tight_layout()
plt.savefig(f'{save_dir}{png_name}_pH{target_pH}{suffix}.png', dpi=300, bbox_inches='tight')
print(f"Saved plots to {png_name}{suffix}.png and {png_name}_pH{target_pH}{suffix}.png")

if args.show:
    plt.show()