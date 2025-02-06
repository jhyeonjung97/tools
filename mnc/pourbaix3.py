import os
import re
import glob
import numpy as np
import pandas as pd
from ase.io import read
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
   
fig_width_pt = 1.8 * 246.0  # LaTeX column width
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # Width in inches
fig_height = fig_width * golden_mean  # Height in inches
fig_size = [fig_width, fig_height]

font_size = 10
tick_font_size = 10

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = 2 * font_size
plt.rcParams['legend.fontsize'] = font_size
plt.rcParams['xtick.labelsize'] = tick_font_size
plt.rcParams['ytick.labelsize'] = tick_font_size
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['lines.linewidth'] = 1.0

# Constants and variables for plotting
Umin, Umax = -0.5, 2.5
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485

pH = np.arange(0, 14, 0.10)
U = np.arange(Umin, Umax, 0.01)
Umax2 = Umax + 0.06 * 14
U2 = np.arange(Umin, Umax2, 0.01)

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.560
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

color = [
    # 'cornflowerblue', ## Sv
    'lightsteelblue', ## S0
    'darkgray', ## S1
    'yellowgreen', # S2
    'red', # S3
    'tan', ## S4
    'pink', ## S5
    'green', # S6
    'salmon', ## S7  
    'blue', # S8
    'khaki', ## S9
    'salmon', ## S10
    'plum', ## S11
]
pH2 = np.arange(0, 14.01, 0.01)

path = '/Users/hailey/Desktop/jan20/pourbaix/'
metal_path = '/Users/hailey/Desktop/jan20/metals.tsv'
# path = '/Users/jiuy97/Desktop/jan6/figures/pourbaix/'
# metal_path = '/Users/jiuy97/Desktop/jan21/metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
metals = ['Fe', 'Co', 'Mo']

elements_data = {
    "Ti": {"electrode_potential": -1.63, "cation_charge": 2},  # Ti^2+ + 2e- → Ti (Ti^2+ prioritized over Ti^4+)
    "V":  {"electrode_potential": -1.18, "cation_charge": 2},  # V^2+ + 2e- → V (V^2+ prioritized over V^3+)
    "Cr": {"electrode_potential": -0.91, "cation_charge": 2},  # Cr^2+ + 2e- → Cr (Cr^2+ prioritized over Cr^3+)
    "Mn": {"electrode_potential": -1.18, "cation_charge": 2},  # Mn^2+ + 2e- → Mn
    "Fe": {"electrode_potential": -0.44, "cation_charge": 2},  # Fe^2+ + 2e- → Fe (Fe^2+ prioritized over Fe^3+)
    "Co": {"electrode_potential": -0.28, "cation_charge": 2},  # Co^2+ + 2e- → Co (Co^2+ prioritized over Co^3+)
    "Ni": {"electrode_potential": -0.25, "cation_charge": 2},  # Ni^2+ + 2e- → Ni
    "Cu": {"electrode_potential": 0.34,  "cation_charge": 2},  # Cu^2+ + 2e- → Cu

    "Zr": {"electrode_potential": -1.45, "cation_charge": 2},  # Zr^2+ + 2e- → Zr (Zr^2+ prioritized over Zr^4+)
    "Nb": {"electrode_potential": -1.00, "cation_charge": 3},  # Nb^3+ + 3e- → Nb (Nb^3+ prioritized over Nb^5+)
    "Mo": {"electrode_potential": -0.20, "cation_charge": 3},  # Mo^3+ + 3e- → Mo (Mo^3+ prioritized over Mo^6+)
    "Tc": {"electrode_potential": 0.74,  "cation_charge": 7},  # Tc^7+ + 7e- → Tc (No lower oxidation state available)
    "Ru": {"electrode_potential": 0.45,  "cation_charge": 2},  # Ru^2+ + 2e- → Ru (Ru^2+ prioritized over Ru^3+)
    "Rh": {"electrode_potential": 0.76,  "cation_charge": 2},  # Rh^2+ + 2e- → Rh (Rh^2+ prioritized over Rh^3+)
    "Pd": {"electrode_potential": 0.95,  "cation_charge": 2},  # Pd^2+ + 2e- → Pd

    "Hf": {"electrode_potential": -1.55, "cation_charge": 4},  # Hf^4+ + 4e- → Hf (No lower oxidation state available)
    "Ta": {"electrode_potential": -1.10, "cation_charge": 3},  # Ta^3+ + 3e- → Ta (Ta^3+ prioritized over Ta^5+)
    "W":  {"electrode_potential": -0.11, "cation_charge": 4},  # W^4+ + 4e- → W (W^4+ prioritized over W^6+)
    "Re": {"electrode_potential": 0.30,  "cation_charge": 4},  # Re^4+ + 4e- → Re (Re^4+ prioritized over Re^7+)
    "Os": {"electrode_potential": 0.40,  "cation_charge": 2},  # Os^2+ + 2e- → Os (Os^2+ prioritized over Os^4+)
    "Ir": {"electrode_potential": 1.16,  "cation_charge": 2},  # Ir^2+ + 2e- → Ir (Ir^2+ prioritized over Ir^3+)
    "Pt": {"electrode_potential": 1.20,  "cation_charge": 2},  # Pt^2+ + 2e- → Pt
}

def addH(x, y):
    return -gh + dgh + 1 * (y + x * const)
    
def addO(x, y):
    return -go + dgo - 2 * (y + x * const)

def addOH(x, y):
    return -goh + dgoh - (y + x * const)

def addOOH(x, y):
    return -gooh + dgooh - 3 * (y + x * const)

def dg(i, x, y):
    if surfs[i][0] is None:
        return None
    elif i == 0 and surfs[i][1] == 2:
        return (surfs[i][0] 
                - surfs[1][0] 
                + bulk_metal
                + surfs[i][1] * addH(x, y) 
                + surfs[i][2] * addO(x, y) 
                + surfs[i][3] * addOH(x, y) 
                + surfs[i][4] * addOOH(x, y))
    return (surfs[i][0] 
            - surfs[1][0] 
            + surfs[i][1] * addH(x, y) 
            + surfs[i][2] * addO(x, y) 
            + surfs[i][3] * addOH(x, y) 
            + surfs[i][4] * addOOH(x, y))
  
for i, metal in enumerate(metals):
    
    A, B = i+1, metal                
    bulk_metal = metal_df.at[B, 'energy']
    
    df = pd.DataFrame()
    df_path = os.path.join(path, f'{A}{B}_energies.tsv')
    df = pd.read_csv(df_path, delimiter='\t', index_col=0)
    df.loc['vac', 'E'] = -.28183609E+03            
    df.loc['vac', ['#H', '#O', '#OH', '#OOH']] = [2, 0, 0, 0]
    
    OER = pd.DataFrame()
    OER_path = os.path.join(path, f'{A}{B}_oer.tsv')
    OER = pd.read_csv(OER_path, delimiter='\t', index_col=0)

    ORR = pd.DataFrame()
    ORR_path = os.path.join(path, f'{A}{B}_orr.tsv')
    ORR = pd.read_csv(ORR_path, delimiter='\t', index_col=0)
    
    surfs = [
        df.loc['vac', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        df.loc['o-oh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),
        # df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
        # df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].tolist(),  
    ]
    
    if A == '3' and B == 'Mo':
        surfs.append(df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].tolist())
    
    surfs = [surf for surf in surfs if not any(pd.isna(x) for x in surf)]
    
    nsurfs = len(surfs)
    lowest_surfaces = []
    
    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs) if dg(k, 0, j) is not None]
        lowest_surfaces.append(np.argmin(values))
        
    crossover = []
    uniquesurf = [lowest_surfaces[0]]
    old_value = lowest_surfaces[0]
    crossover.append(Umin)
    
    for j in range(len(U2)):
        if lowest_surfaces[j] != old_value:
            uniquesurf.append(lowest_surfaces[j])
            crossover.append(U2[j])
            old_value = lowest_surfaces[j]
    crossover.append(Umax2)
    
    # plt.clf()
    # fig = plt.figure(figsize=fig_size, dpi=300)
    # ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    # ax.axis([0, 14, Umin, Umax])
    # ax.set_xlabel(r'pH', fontsize='large')
    # ax.set_ylabel(r'U/V', fontsize='large')
    # current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    # extraticks = [1.23]
    # combined_ticks = sorted(set(current_yticks) | set(extraticks))
    # plt.yticks(combined_ticks)
    # plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # for i in range(len(uniquesurf)):
    #     k = uniquesurf[i]
    #     label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
    #     plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
    #                      facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
    #     plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)
    # plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    # if A == 1 and B == 'Fe':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][1] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.70, 
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][1] - 0.96, # 1.14, 
    #             r"S$_4$$\rightarrow$S$_7$$\rightarrow$S$_{10}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][1]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='gray', rotation=-9.5, fontsize=10)
    # elif A == 2 and B == 'Co':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.70,
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][2] - 0.98, 
    #             r"S$_5$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][2]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='gray', rotation=-9.5, fontsize=10)
    # elif A == 3 and B == 'Mo':
    #     plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][2] - pH2 * const, '--', color='green', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.8, OER['onsetP'][2] - 1.30, 
    #             r"S$_5$$\rightarrow$S$_8$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$: " + f"{OER['overP'][2]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][2] - 0.34,
    #             r"S$_{11}$$\rightarrow$S$_{10}$$\rightarrow$S$_8$$\rightarrow$S$_5$: " + f"{ORR['overP'][2]:.2f} eV", 
    #             color='green', rotation=-9.5, fontsize=10)
    # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
    #            ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
    #            fancybox=True, shadow=True)
    # plt.savefig(f'{path}{A}{B}_pourbaix_full.png', bbox_inches='tight')
    # print(f"Figure saved as {A}{B}_pourbaix_full.png")
    # plt.close()

    # plt.clf()
    # fig = plt.figure(figsize=fig_size, dpi=300)
    # ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    # ax.axis([0, 14, Umin, Umax])
    # ax.set_xlabel(r'pH', fontsize='large')
    # ax.set_ylabel(r'U/V', fontsize='large')
    # current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    # extraticks = [1.23]
    # combined_ticks = sorted(set(current_yticks) | set(extraticks))
    # plt.yticks(combined_ticks)
    # plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # for i in range(len(uniquesurf)):
    #     k = uniquesurf[i]
    #     plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
    #                      facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
    #     plt.plot([], [], color=color[k], alpha=0.3, linewidth=5)
    # plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    # if A == 1 and B == 'Fe':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][1] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.70, 
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][1] - 0.96, # 1.14, 
    #             r"S$_4$$\rightarrow$S$_7$$\rightarrow$S$_{10}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][1]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='gray', rotation=-9.5, fontsize=10)
    # elif A == 2 and B == 'Co':
    #     plt.plot(pH2, OER['onsetP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][0] - pH2 * const, '--', color='gray', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][0] - 0.70,
    #             r"S$_1$$\rightarrow$S$_4$$\rightarrow$S$_5$$\rightarrow$S$_8$: " + f"{OER['overP'][0]:.2f} eV", 
    #             color='black', rotation=-9.5, fontsize=10)
    #     ax.text(6.5, OER['onsetP'][2] - 0.98, 
    #             r"S$_5$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$$\rightarrow$S$_{13}$: " + f"{OER['overP'][2]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][0] - 0.58, 
    #             r"S$_8$$\rightarrow$S$_5$$\rightarrow$S$_4$$\rightarrow$S$_1$: " + f"{ORR['overP'][0]:.2f} eV", 
    #             color='gray', rotation=-9.5, fontsize=10)
    # elif A == 3 and B == 'Mo':
    #     plt.plot(pH2, OER['onsetP'][2] - pH2 * const, '--', color='red', lw=1, dashes=(3, 1))
    #     plt.plot(pH2, ORR['onsetP'][2] - pH2 * const, '--', color='green', lw=1, dashes=(3, 1))
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(6.8, OER['onsetP'][2] - 1.30, 
    #             r"S$_5$$\rightarrow$S$_8$$\rightarrow$S$_{10}$$\rightarrow$S$_{11}$: " + f"{OER['overP'][2]:.2f} eV", 
    #             color='red', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
    #     ax.text(0.2, ORR['onsetP'][2] - 0.34,
    #             r"S$_{11}$$\rightarrow$S$_{10}$$\rightarrow$S$_{8}$$\rightarrow$S$_{5}$: " + f"{ORR['overP'][2]:.2f} eV", 
    #             color='green', rotation=-9.5, fontsize=10)
    # plt.savefig(f'{path}{A}{B}_pourbaix_clean.png', bbox_inches='tight')
    # print(f"Figure saved as {A}{B}_pourbaix_clean.png")
    # # plt.show()
    # plt.close()
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([-1.0, 2.5, -600, 200])
    ax.set_xlabel(r'RHE (V)', fontsize='large')
    ax.set_ylabel(r'$\Delta$G (kJ/mol)', fontsize='large')
    xx = np.arange(-1.00, 2.55, 0.01)
    for k in range(nsurfs):
        label = r"S$_{%i}$(H: %i O: %i OH: %i OOH: %i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        dg_value = dg(k, 0, xx)
        if dg_value is not None:
            ax.plot(xx, dg_value * kjmol, '-', lw=1, c=color[k], label=label)
        else:
            print(f"Skipping plot for surface {k} due to missing data.")
    plt.xlim(-1.0, 2.5)
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    plt.savefig(f'{path}{A}{B}_pourbaix.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix.png")
    # plt.show()
    plt.close()
    