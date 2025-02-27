import os
from ase.io import read, write
import numpy as np # v1.26.4
import pandas as pd # v2.2.2
import matplotlib.pyplot as plt # v3.8.3
from matplotlib.lines import Line2D

# Scaling relationship
a1, b1 = 1.79, 0.62  # OH → O
a2, b2 = 0.87, 3.08  # OH → OOH

# Gas phase energies
h2 = -6.77149190
h2o = -14.23091949

zpeh2o, cvh2o, tsh2o = 0.558, 0.103, 0.675
zpeh2, cvh2, tsh2 = 0.268, 0.0905, 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
go = gh2o - gh2 - 1.229 * 2

# Adsorbate corrections
zpeoh, cvoh, tsoh = 0.376, 0.042, 0.066
zpeo, cvo, tso = 0.064, 0.034, 0.060
zpeooh, cvooh, tsooh = 0.471, 0.077, 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh

metals = {'Cr': 'blue', 'Mn': 'orange', 'Fe': 'green', 'Co': 'red', 'Ni': 'purple'}
adsorbates = ['O', 'OH']

def main():
    df = pd.DataFrame(columns=['a.site', 'mCr', 'mMn', 'mFe', 'mCo', 'mNi'])    
    for i, (metal, color) in enumerate(metals.items()):
        for j in range(3):
            path = f'/pscratch/sd/j/jiuy97/5_HEO/4_local/{i+1}_{metal}/{j+1}_'
            path_DONE = os.path.join(path, 'DONE')
            path_json = os.path.join(path, 'final_with_calculator.json')
            if os.path.exists(path_DONE) and os.path.exists(path_json):
                df.at[f'{metal}{j+1}', 'a.site'] = metal
                atoms = read(path_json)
                moments = atoms.get_magnetic_moments()
                for m in metals:
                    for atom in atoms:
                        if atom.symbol == m and not atom.symbol == metal:   
                            df.at[f'{metal}{j+1}', f'm{m}'] = moments[atom.index]
                df.at[f'{metal}{j+1}', f'm{metal}'] = moments[17]
                df.at[f'{metal}{j+1}', 'clean'] = atoms.get_total_energy()
            else:
                df.at[f'{metal}{j+1}', 'clean'] = np.nan
                
            for k, ads in enumerate(adsorbates):
                path = f'/pscratch/sd/j/jiuy97/5_HEO/4_local/{i+1}_{metal}/{j+1}_/{k+1}_{ads}'
                path_DONE = os.path.join(path, 'DONE')
                path_json = os.path.join(path, 'final_with_calculator.json')
                if os.path.exists(path_DONE) and os.path.exists(path_json):
                    atoms = read(path_json)
                    df.at[f'{metal}{j+1}', ads] = atoms.get_total_energy()
                else:
                    df.at[f'{metal}{j+1}', ads] = np.nan
                    
    for index, row in df.iterrows():
        if pd.notna(row["clean"]) and pd.notna(row["O"]):
            df.at[index, "go"] = row["O"] - row["clean"] - (gh2o - gh2) + dgo
        if pd.notna(row["clean"]) and pd.notna(row["OH"]):
            df.at[index, "goh"] = row["OH"] - row["clean"] - (gh2o - gh2/2) + dgoh
            
    df = df.apply(pd.to_numeric, errors='ignore')
    df.to_csv(f'/pscratch/sd/j/jiuy97/5_HEO/4_local/scaling.csv', sep=',')
    df.to_csv(f'/pscratch/sd/j/jiuy97/5_HEO/4_local/scaling.tsv', sep='\t', float_format='%.2f')

    df['color'] = df['a.site'].map(metals)
    scaling_OH_O(df)

def scaling_OH_O(df):
    """Plots the ORR volcano plot."""
    xmin, xmax, xtick = -1.0, 4.0, 0.5
    ymin, ymax, ytick = 0.0, 5.0, 0.5

    x, y, c = df['goh'], df['go'], df['color']
    
    plt.figure(figsize=(4, 3), dpi=200)
    plt.scatter(x, y, c=c, s=20, zorder=3)
    
    xx = np.linspace(min(x), max(x), 100)
    coeffs = np.polyfit(x, y, 1)
    line = np.poly1d(coeffs)
    plt.plot(xx, line(xx), linewidth=0.5, linestyle='--', color='black', zorder=2)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xticks(np.arange(xmin, xmax+xtick, xtick), fontsize=8) 
    plt.yticks(np.arange(ymin, ymax+ytick, ytick), fontsize=8) 
    plt.xlabel(r"$\mathrm{\Delta G_{OH}}$ (eV)")
    plt.ylabel(r"$\mathrm{\Delta G_{O}}$ (eV)")

    legend_elements = [
        Line2D([0], [0], marker='o', linestyle='', label='Cr', color='w', markerfacecolor='blue'),
        Line2D([0], [0], marker='o', linestyle='', label='Mn', color='w', markerfacecolor='orange'),
        Line2D([0], [0], marker='o', linestyle='', label='Fe', color='w', markerfacecolor='green'),
        Line2D([0], [0], marker='o', linestyle='', label='Co', color='w', markerfacecolor='red'),
        Line2D([0], [0], marker='o', linestyle='', label='Ni', color='w', markerfacecolor='purple'),
    ]
    plt.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, 0.02), 
               fontsize=8, ncol=5, columnspacing=0.5, handletextpad=0.0)
    plt.tight_layout()
    plt.savefig('heo_scaling.png')
    plt.show()

if __name__ == "__main__":
    main()
