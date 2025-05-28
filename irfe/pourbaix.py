import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase import Atoms
import matplotlib.colors as mcolors
import pandas as pd

# 전압 범위 정의
Umin, Umax = -0.5, 2.0
ymin, ymax = -1.0, 1.0

def count_elements(atoms):
    """원자 구조에서 O와 H의 개수를 세는 함수"""
    symbols = atoms.get_chemical_symbols()
    return {'O': symbols.count('O'), 'H': symbols.count('H')}

def read_vib_correction(vib_file):
    """vib.txt 파일에서 G(T) 보정값을 읽는 함수"""
    try:
        with open(vib_file, 'r') as f:
            content = f.read()
            for line in content.split('\n'):
                if 'Zero-point energy E_ZPE' in line:
                    zpe = float(line.split()[-2])
                if 'Entropy contribution T*S' in line:
                    ts = float(line.split()[-2])
            return zpe - ts
    except:
        return 0.0

def calculate_free_energy(E, E0, n_H, n_O, U, G_vib=0.0, pH=0):
    """Gibbs free energy 계산 - 진동 자유에너지와 열역학적 보정 포함"""
    kT = 0.026
    const = kT * np.log(10)
    
    e_h2 = -6.77108058
    e_h2o = -14.22266334
    
    zpe_h2o = 0.558
    cv_h2o = 0.103
    ts_h2o = 0.675

    zpe_h2 = 0.268
    cv_h2 = 0.0905
    ts_h2 = 0.408

    g_h2o = e_h2o + zpe_h2o - ts_h2o + cv_h2o
    g_h2 = e_h2 + zpe_h2 - ts_h2 + cv_h2

    surface_term = E - E0 + G_vib
    n_e = (n_H - 2*n_O)
    
    G = surface_term + n_e * U + n_e * const * pH
    G -= n_H * 0.5 * g_h2
    G -= n_O * (g_h2o - g_h2)
    
    return G/8

def create_pourbaix_diagram(base_path):
    metal_systems = {
        '0_Ir': 'Ir',
        '1_Mn': 'Mn',
        '2_Fe': 'Fe',
        '3_Co': 'Co',
        '4_Ni': 'Ni'
    }

    coverage_settings = {
        '1_H': '2_layer_hol',
        '2_OH': '2_layer_brg',
        '3_O': '2_layer_hol'
    }

    color_dict = {
        'clean': 'black',
        'H': 'blue',
        'OH': 'green',
        'O': 'red'
    }

    U_range = [Umin, Umax]

    for m, target_metal in enumerate(metal_systems.keys()):
        clean_path = os.path.join(base_path, 'slab', target_metal, 'final_with_calculator.json')
        if not os.path.exists(clean_path):
            continue
        clean_done_path = os.path.join(base_path, 'slab', target_metal, 'DONE')
        if not os.path.exists(clean_done_path):
            continue

        clean_atoms = read(clean_path)
        energy_clean = clean_atoms.get_potential_energy()

        stable_sites = {}

        for ads, coverage in coverage_settings.items():
            metal_path = os.path.join(base_path, ads, target_metal)
            if not os.path.exists(metal_path):
                continue

            if target_metal == '0_Ir' and ads == '1_H':
                site_path = os.path.join(metal_path, '1_layer_top')
            else:
                site_path = os.path.join(metal_path, coverage)
            # site_path = os.path.join(metal_path, coverage)

            if not os.path.exists(site_path):
                continue

            full_path = os.path.join(site_path, 'final_with_calculator.json')
            if not os.path.exists(full_path):
                continue

            done_path = os.path.join(site_path, 'DONE')
            if not os.path.exists(done_path):
                continue

            vib_file = os.path.join(site_path, 'vib.txt')
            G_vib = read_vib_correction(vib_file)
            atoms = read(full_path)
            energy = atoms.get_potential_energy()
            counts = count_elements(atoms)

            Gmin = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umin, G_vib)
            Gmax = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umax, G_vib)
            G0 = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], 0, G_vib)

            ads_type = ads.split('_')[1]
            stable_sites[ads_type] = {
                'Gmin': Gmin,
                'Gmax': Gmax,
                'G0': G0
            }

        plt.figure(figsize=(5, 4))
        plt.plot(U_range, [0]*len(U_range), color='black', label='clean')
        plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

        print(f"\n{metal_systems[target_metal]}:")
        
        for ads_type, entry in stable_sites.items():
            color = color_dict.get(ads_type, 'gray')
            plt.plot(U_range, [entry['Gmin'], entry['Gmax']], color=color, label=ads_type)
            
            # 기울기와 y절편 계산
            slope = (entry['Gmax'] - entry['Gmin']) / (Umax - Umin)
            y_intercept = entry['Gmin'] - slope * Umin
            
            # x절편 계산 (y = 0일 때의 x값)
            x_intercept = -y_intercept / slope if slope != 0 else None
            
            print(f"{ads_type}\t{x_intercept:.2f} V\t{y_intercept:.2f} eV")

        plt.xlabel('Potential (V vs RHE)')
        plt.ylabel('Relative Gibbs Free Energy (eV)')
        plt.xlim(Umin, Umax)
        plt.ylim(ymin, ymax)
        plt.legend(bbox_to_anchor=(1.05, 1.05), loc='upper left', fontsize='small', labelspacing=0.2)
        
        save_path = os.path.join(base_path, 'figures', f'pourbaix_diagram_{m}{metal_systems[target_metal]}.png')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    base_path = "/Users/hailey/Desktop/4_IrFe3"
    create_pourbaix_diagram(base_path) 