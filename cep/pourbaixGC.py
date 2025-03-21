### Hailey ###
### Date: 2025-03-13 ###

import os
import numpy as np
from math import log10
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

GCDFT = True
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
    # [-280.1784, 1, +0, 0, 0, 0, 0, -0.3714,  0.2173, -1.96, 'clean'],
    # [-282.5588, 1, +0, 1, 0, 0, 0, -0.3438, -0.3362, -0.89, '*H(Fe)'],
    # [-282.6377, 1, +0, 1, 0, 0, 0, -0.3608, -0.6948, -0.83, '*H(N)'],
    # [-290.9780, 1, +0, 0, 1, 0, 0, -0.8300, -0.4206, -4.45, '*OH'],
    # [-285.9474, 1, +0, 0, 0, 1, 0, -0.7036,  0.3162, -5.73, '*O'],
    # [-300.7382, 1, +0, 0, 2, 0, 0, -0.7214,  0.3152, -5.60, '*OH+*OH(adg)'],
    # [-300.6264, 1, +0, 0, 2, 0, 0, -0.9681, -0.9980, -3.54, '*OH+*OH(anti)'],
    # [-295.2283, 1, +0, 0, 1, 1, 0, -0.9971, -1.1225, -3.47, '*OH+*O(adg)'],
    # [-295.5266, 1, +0, 0, 1, 1, 0, -0.8134,  0.0605, -4.80, '*OH+*O(anti)'],
    # [-289.6763, 1, +0, 0, 0, 2, 0, -0.9297, -0.7491, -3.28, '*O+*O(anti)'],
]

nsurfs = len(surfs)
surface_reference = surfs[0][0]
for k in range(nsurfs):
    surfs[k][9] += surfs[k][0] ## remove this if you have GCDFT data
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

# Set Axes Limits and Labels
fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
ax.axis([0, 14, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

# 색상 정의
colors = list(colormaps["tab20c"].colors) + list(colormaps["tab20b"].colors)

# lowest_surfaces가 실제로 어떤 값들로 되어 있는지 확인하고, 고유값 리스트 생성
unique_ids = np.unique(lowest_surfaces)
nsurfs = len(unique_ids)

# 색상 수보다 surface 수가 많으면 에러 방지
if nsurfs > len(colors):
    raise ValueError("Surface 종류가 너무 많아 색상 부족!")

# 고유 surface ID를 0부터 차례로 다시 매핑 (ex: {10:0, 15:1, 30:2, ...})
id_map = {val: idx for idx, val in enumerate(unique_ids)}
mapped_surfaces = np.vectorize(id_map.get)(lowest_surfaces)

# 컬러맵과 정규화 설정
cmap = mcolors.ListedColormap(colors[:nsurfs])
bounds = np.arange(nsurfs + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# legend 생성 - 고유 surface ID 순서에 따라
print(unique_ids)
for idx, surf_id in enumerate(unique_ids):
    label = surfs[int(surf_id)][10]
    plt.plot([], [], color=colors[idx], linewidth=5, label=label)

# pcolormesh
pH, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH, U, mapped_surfaces, cmap=cmap, norm=norm)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1,
       fontsize='x-small', handlelength=3, edgecolor='black')

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