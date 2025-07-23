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
import pandas as pd
from ase.db import connect
import re
from collections import Counter

# Cobalt-specific constants and energies (from pourbaix-co-wide.py)
bulk_metal = -3.74  # Co, eV (gc in pourbaix-co-wide.py, per atom)
concentration_ion = 1e-6
concentration_gas = 1e-6

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram for Cobalt')
parser.add_argument('--gc', action='store_true', help='Enable GCDFT mode')
parser.add_argument('--bulk', action='store_true', help='Enable bulk Pourbaix mode')
parser.add_argument('--suffix', type=str, default='CoSAC', help='Suffix for output filename')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--save-dir', action='store_true', help='Save to predefined directory')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--x', type=float, default=6, help='Figure width in inches (default: 6)')
parser.add_argument('--y', type=float, default=7, help='Figure height in inches (default: 7)')
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

# units
kjmol = 96.485
calmol = 23.061

# constants
kb = 8.617e-5 # eV/K
T = 298.15 # K
const = kb * T * np.log(10) # 0.0592 eV
water = 2.4583 # eV (from pourbaix-co-wide.py)

# ticks
tick = args.tick
pHmin, pHmax = -2, 16
pHrange = np.arange(pHmin, pHmax + tick, tick)
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

gn2 = -16.64503942 + 0.098 - 0.592  # n2 + zpen2 - tsn2

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

# Cobalt ions (from get_ion_entries in pourbaix-co-wide.py)
ions = [
    # ['Ef', '#M(Co)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'concentration', 'name']
    [ -12.800/calmol, 1, +2, 0, 0, 0, 0, 0, 0, 0, concentration_ion, 'Co²⁺(aq)'],
    [  28.900/calmol, 1, +3, 0, 0, 0, 0, 0, 0, 0, concentration_ion, 'Co³⁺(aq)'],
    [ -82.970/calmol, 1, +1, 0, 2, 0, 0, 0, 0, concentration_ion, 'HCoO₂⁻(aq)'],
]

# Cobalt solids (from get_solid_entries in pourbaix-co-wide.py)
solids = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 'Co(s)'],
    [-52.310/calmol, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 'CoO(s)'],
    [-167.835/calmol, 3, 0, 0, 0, 4, 0, 0, 0, 0, 1, 'Co₃O₄(s)'],
    [-51.840/calmol, 1, 0, 0, 0, 2, 0, 0, 0, 0, 1, 'CoO₂(s)'],
    [-109.999/calmol, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 'Co(OH)₂(s)'],
    [-142.600/calmol, 1, 0, 0, 3, 0, 0, 0, 0, 0, 1, 'Co(OH)₃(s)'],
]

def parse_formula(formula):
    # 예: 'CoC26H2N4O2' -> {'Co':1, 'C':26, 'H':2, 'N':4, 'O':2}
    pat = re.compile(r'([A-Z][a-z]*)([0-9]*)')
    d = {}
    for el, num in pat.findall(formula):
        d[el] = int(num) if num else 1
    return d

db = connect('combined.db')
rows = list(db.select())
surfs = []

for i, row1 in enumerate(rows):
    f1 = parse_formula(row1.formula)
    for j, row2 in enumerate(rows):
        if i == j:
            continue
        f2 = parse_formula(row2.formula)
        # 공통 부분
        common = {el: min(f1.get(el,0), f2.get(el,0)) for el in set(f1) | set(f2)}
        # 유니크 부분 (row1 기준)
        unique = {el: f1[el] - common.get(el,0) for el in f1 if f1[el] - common.get(el,0) > 0}
        # 예시: CoC26H2N4O2 vs CoC26N4O2 -> unique = {'H':2}
        # 필요한 원소만 추출해서 surfs에 넣기
        nCo = unique.get('Co', 0)
        nN = unique.get('N', 0)
        nH = unique.get('H', 0)
        nO = unique.get('O', 0)
        # ... 필요에 따라 OH, OOH 등 추가
        surfs.append([row1.energy, nCo, nN, 0, nH, 0, nO, 0, 0, 0, 1, row1.formula + '-' + row2.formula])

print(surfs)
nsurfs = len(surfs)
ref_surf = surfs[0][0]
ref_surf_gc = surfs[0][0]
for k in range(nsurfs):
    formation_energy_corr = (
        - surfs[k][3] * (gh2 - dgh) # H
        - surfs[k][4] * (goh - dgoh) # OH
        - surfs[k][5] * (go - dgo) # O
        - surfs[k][6] * (gooh - dgooh) # OOH
    )
    surfs[k][0] = surfs[k][0] - ref_surf + formation_energy_corr 
    surfs[k][10] = surfs[k][0]
    surfs[k][11] = 1

def dg(k, pH, U, n_ref):
    if GCDFT:
        surface_term = ((surfs[k][8]*(U**2) + surfs[k][9]*U + surfs[k][10])
                          - (surfs[n_ref][8]*(U**2) + surfs[n_ref][9]*U + surfs[n_ref][10]))
    else:
        surface_term = surfs[k][0] - surfs[n_ref][0]
    U_coeff = 1*surfs[k][3] - 1*surfs[k][4] - 2*surfs[k][5] - 3*surfs[k][6] - surfs[k][2]
    pH_coeff = 1*surfs[k][3] - 1*surfs[k][4] - 2*surfs[k][5] - 3*surfs[k][6]
    dg = surface_term + U_coeff*U + const*pH_coeff*pH + const*log10(surfs[k][11])
    return dg

# The rest of the code (plotting, finding lowest surfaces, etc.) follows the same structure as ti-sac.py
# ... existing code ... 