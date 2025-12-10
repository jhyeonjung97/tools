import os
import pandas as pd
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["axes.prop_cycle"] = plt.cycler(color=plt.get_cmap("tab20").colors)

# gas (전역 상수: 두 스크립트에서 동일 값 사용)
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
gh = gh2/2
go = gh2o - gh2
goh = gh2o - gh
goo = 2*gh2o - 2*gh2
gooh = 2*gh2o - 3*gh

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
dgoo = dgooh - dgh

def extract_energy_from_json(json_path):
    """
    JSON 파일에서 ASE Atoms 객체를 읽고 에너지를 추출하는 함수
    
    Args:
        json_path (str): JSON 파일 경로
        
    Returns:
        float: 에너지 값 (eV)
    """
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def main():
    plt.figure(figsize=(4, 4))
    plt.plot([0, 2], [0, 2], color='black', linewidth=0.5, zorder=2)
    plt.plot([0, 2], [3.2, 1.2], color='black', linewidth=0.5, zorder=2)         
    plt.plot([0, 2], [2.46, 0.46], color='black', linewidth=0.5, zorder=2)
    plt.axhline(1.23, color='black', linewidth=0.5, zorder=1, dashes=(5, 5))
    plt.axvspan(1.23, 1.60, color='gainsboro', zorder=0)    # mistyrose
    plt.xlabel(r'$\Delta G_{O} - \Delta G_{OH}$ (eV)')
    plt.ylabel(r'$\Delta G_{max}$ (eV)')
    plt.xlim(0.6, 2)
    plt.ylim(1.1, 2.1)
    plt.grid(True, alpha=0.1)
    plt.gca().invert_yaxis()
    plt.savefig("volcano-scheme.png", dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
