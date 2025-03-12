import os
import numpy as np
from math import log10
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

GCDFT = False
bulk_metal = -5.7954 # Fe, eV

# constants
kb = 8.617e-5 # eV/K
T = 298.15 # K
const = kb * T * np.log(10) # 0.0592 eV

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
go = gh2o - gh2 - 2.46
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
    [-55.910/calmol,  1, +2, 1, 0, 1, 0, 0, 0, 0, 'FeOH²⁺(aq)'],
    [-106.200/calmol, 1, +1, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂⁺(aq)'],
]

solids = [
    # ['Ef', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [0,               1, +0, 0, 0, 0, 0, 0, 0, 0, 'Fe(s)'],
    [-58.880/calmol,  1, +0, 0, 0, 1, 0, 0, 0, 0, 'FeO'],
    [-242.400/calmol, 3, +0, 0, 0, 4, 0, 0, 0, 0, 'Fe₃O₄'],
    [-177.100/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-161.930/calmol, 2, +0, 0, 0, 3, 0, 0, 0, 0, 'Fe₂O₃'],
    [-115.570/calmol, 1, +0, 0, 2, 0, 0, 0, 0, 0, 'Fe(OH)₂'],
    [-166.000/calmol, 1, +0, 0, 3, 0, 0, 0, 0, 0, 'Fe(OH)₃'],
]

nions, nsolids = len(ions), len(solids)
for i in range(nions):
    if ions[i][1] > 1:
        ions[i] = [x / ions[i][1] if isinstance(x, (int, float)) else x for x in ions[i]]
for s in range(nsolids):
    if solids[s][1] > 1:
        solids[s] = [x / solids[s][1] if isinstance(x, (int, float)) else x for x in solids[s]]
    
surfs = [
    # ['E', '#M(=Fe)', '#e', '#H', '#OH', '#O', '#OOH', 'A', 'B', 'C', 'name']
    [-271.9532, 0, +0, 0, 0, 0, 0, -0.3442, -0.1279, 0, 'vac'],
    [-281.8361, 0, +0, 2, 0, 0, 0, -0.3342, -0.1079, 0, 'vac(H₂)'],
    [-280.1784, 1, +0, 0, 0, 0, 0, -0.3660,  0.0665, 0, 'clean'],
    [-282.5588, 1, +0, 1, 0, 0, 0, -0.3714,  0.2173, 0, '＊H(Fe)'],
    [-282.6377, 1, +0, 1, 0, 0, 0, -0.3438, -0.3362, 0, '＊H(N)'],
    [-290.9780, 1, +0, 0, 1, 0, 0, -0.8300, -0.4206, 0, '＊OH'],
    [-285.9474, 1, +0, 0, 0, 1, 0, -0.3608, -0.6948, 0, '＊O'],
    [-300.7382, 1, +0, 0, 2, 0, 0, -0.7036,  0.3162, 0, '＊OH+＊OH(adg)'],
    [-300.6264, 1, +0, 0, 2, 0, 0, -0.7214,  0.3152, 0, '＊OH+＊OH(anti)'],
    [-295.2283, 1, +0, 0, 1, 1, 0, -0.8385, -0.6094, 0, '＊OH+＊O(adg)'],
    [-295.5266, 1, +0, 0, 1, 1, 0, -0.7100,  0.1638, 0, '＊OH+＊O(anti)'],
    [-289.6763, 1, +0, 0, 0, 2, 0, -0.6429,  2.1679, 0, '＊O+＊O(anti)'],
]
surfs[9] = surfs[0] ## remove this if you have GCDFT data

nsurfs = len(surfs)
for k in range(nsurfs):
    formation_energy_corr = (
        - surfs[0][0] # surface_ref(vac)
        - surfs[k][3] * (gh - dgh) # H
        - surfs[k][4] * (go - dgo) # O
        - surfs[k][5] * (goh - dgoh) # OH
        - surfs[k][6] * (gooh - dgooh) # OOH
    )
    surfs[k][0] += formation_energy_corr # E
    surfs[k][9] += formation_energy_corr # C
            
new_surfs = []
for k in range(nsurfs):
    if surfs[k][1] == 0:
        for i in range(nions):
            new_surfs.append(surfs[k] + ions[i])
        for s in range(nsolids):
            new_surfs.append(surfs[k] + solids[s])

surfs.extend(new_surfs)  # Add new surfaces after looping
nsurfs = len(surfs)  # Update length

lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

pHindex = 0
for pH in pHrange:
    Uindex = 0
    for U in Urange:
        values = []
        for k in range(nsurfs):
            values.append(dg(k, pH, U, concentration=1e-6))
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces[Uindex][pHindex] = sorted_values[0]
        Uindex+=1
    pHindex+=1

# Set Axes Limits and Labels
fig = plt.figure(figsize=(8, 6), dpi=100)
ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])
ax.axis([0, 14, Umin, Umax])
ax.set_xlabel('pH', labelpad=0)
ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
ax.tick_params(right=True, direction="in")

colors = plt.cm.tab20.colors[:nsurfs]
cmap = mcolors.ListedColormap(colors)
bounds = np.arange(nsurfs + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

for k in range(nsurfs):
    label = surfs[k][10]
    # label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][3], surfs[k][4], surfs[k][5], surfs[k][6])
    plt.plot([], [], color=colors[k], linewidth=5, label=label)

pH, U = np.meshgrid(pHrange, Urange)
plt.pcolormesh(pH, U, lowest_surfaces, cmap=cmap, norm=norm)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1,
       fontsize='x-small', handlelength=3, edgecolor='black')

plt.plot(pHrange, 1.23-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 1.00, r'2H$_2$O $\leftrightarrow$ 4H$^+$+O$_2$+4e$^-$',
        color='blue', rotation=-8, fontsize=10)
plt.plot(pHrange, 0-pHrange*const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, -0.15 , r'H$_2 $ $\leftrightarrow$ 2H$^+$+$\ $2e$^-$',
        color='blue', rotation=-8, fontsize=10)

fig.savefig(f'pourbaixGC.png', bbox_inches='tight') 