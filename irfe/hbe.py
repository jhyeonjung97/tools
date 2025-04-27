import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
import pandas as pd
import re

# 전압 범위 정의
Umin, Umax = -0.5, 1.5
ymin, ymax = -2.0, 1.0

def count_hydrogen(atoms):
    symbols = atoms.get_chemical_symbols()
    return symbols.count('H')

def calculate_free_energy(E, E0, n_H, U, G_vib=0.0, pH=0):
    kT = 0.026
    const = kT * np.log(10)
    e_h2 = -6.77108058
    zpe_h2 = 0.268
    cv_h2 = 0.0905
    ts_h2 = 0.408
    g_h2 = e_h2 + zpe_h2 - ts_h2 + cv_h2
    surface_term = E - E0 + G_vib
    n_e = n_H
    G = surface_term + n_e * U + n_e * const * pH
    G -= n_H * 0.5 * g_h2
    return G/8

def get_label_from_subfolder(subfolder):
    return subfolder

def get_gvib_dict(csv_path):
    """csv에서 H/layer/top, H/layer/hol의 G_vib 값을 dict로 반환"""
    df = pd.read_csv(csv_path)
    gvib_dict = {}
    for _, row in df.iterrows():
        if row['ads'] == 'H' and row['site1'] == 'layer':
            if row['site2'] == 'top':
                gvib_dict['top'] = row['G_vib']
            elif row['site2'] == 'hol':
                gvib_dict['hol'] = row['G_vib']
    return gvib_dict

def save_data_to_files(data, U_range, base_path, folder, h_folder):
    """데이터를 CSV와 TSV 파일로 저장"""
    df_data = []
    for entry in data:
        Gmin = calculate_free_energy(entry['energy'], entry['energy0'], entry['n_H'], Umin, entry['G_vib'])
        Gmax = calculate_free_energy(entry['energy'], entry['energy0'], entry['n_H'], Umax, entry['G_vib'])
        G0 = calculate_free_energy(entry['energy'], entry['energy0'], entry['n_H'], 0, entry['G_vib'])
        row_data = {
            'label': entry['label'],
            'n_H': entry['n_H'],
            'E': entry['energy'],
            'E0': entry['energy0'],
            'G_vib': entry['G_vib'],
            'Gmin': Gmin,
            'Gmax': Gmax,
            'G0': G0
        }
        df_data.append(row_data)
    df = pd.DataFrame(df_data)
    csv_path = os.path.join(base_path, f"hbe_{folder}_{h_folder}.csv")
    tsv_path = os.path.join(base_path, f"hbe_{folder}_{h_folder}.tsv")
    df.to_csv(csv_path, index=False)
    df.round(2).to_csv(tsv_path, sep='\t', index=False)

def strip_leading_number(s):
    """앞에 붙은 숫자와 언더바를 제거"""
    return re.sub(r'^\d+_', '', s)

def calc_h_adsorption_energies(data, g_h2):
    """
    data: nH별로 정렬된 리스트, 각 원소는 dict (label, n_H, G0 등 포함)
    g_h2: H2의 깁스 자유에너지
    """
    # nH 오름차순 정렬
    sorted_data = sorted(data, key=lambda x: int(re.sub(r'\D', '', x['label'])))
    results = []
    for i in range(1, len(sorted_data)):
        n = sorted_data[i]['n_H']
        G_n = sorted_data[i]['energy'] + sorted_data[i]['G_vib']
        G_n_1 = sorted_data[i-1]['energy'] + sorted_data[i-1]['G_vib']
        delta_G_H = G_n - G_n_1 - 0.5 * g_h2
        results.append({
            'from': sorted_data[i-1]['label'],
            'to': sorted_data[i]['label'],
            'n_H': n,
            'delta_G_H': delta_G_H
        })
    return results

def main():
    base_path = "/Users/hailey/Desktop/4_IrFe3"
    U_range = [Umin, Umax]

    # G_vib 값 미리 읽기
    gvib_ir = get_gvib_dict(os.path.join(base_path, "pourbaix_data_Ir.csv"))
    gvib_irfe = get_gvib_dict(os.path.join(base_path, "pourbaix_data_IrFe.csv"))

    # H2의 깁스 자유에너지 (pourbaix.py와 동일)
    e_h2 = -6.77108058
    zpe_h2 = 0.268
    cv_h2 = 0.0905
    ts_h2 = 0.408
    g_h2 = e_h2 + zpe_h2 - ts_h2 + cv_h2

    for folder in sorted(os.listdir(base_path)):
        if not ("Ir" in folder and os.path.isdir(os.path.join(base_path, folder))):
            continue
        folder_path = os.path.join(base_path, folder)
        folder_stripped = strip_leading_number(folder)

        # 시스템 종류에 따라 clean surface 경로 결정
        if "IrFe" in folder:
            slab_subfolder = "2_Fe"
            gvib_source = gvib_irfe
        else:
            slab_subfolder = "0_Ir"
            gvib_source = gvib_ir

        slab_clean_path = os.path.join(base_path, "slab", slab_subfolder, "final_with_calculator.json")
        if not os.path.exists(slab_clean_path):
            print(f"Clean surface not found: {slab_clean_path}")
            continue
        clean_atoms = read(slab_clean_path)
        energy0 = clean_atoms.get_potential_energy()

        for h_folder in sorted(os.listdir(folder_path)):
            if not ("_H_" in h_folder and os.path.isdir(os.path.join(folder_path, h_folder))):
                continue
            h_folder_path = os.path.join(folder_path, h_folder)
            h_folder_stripped = strip_leading_number(h_folder)

            data = []
            for subfolder in sorted(os.listdir(h_folder_path), key=lambda x: int(re.sub(r'\D', '', x)) if re.match(r'\d+H', x) else 0):
                subfolder_path = os.path.join(h_folder_path, subfolder)
                json_path = os.path.join(subfolder_path, "final_with_calculator.json")
                if not os.path.exists(json_path):
                    continue
                atoms = read(json_path)
                energy = atoms.get_potential_energy()
                n_H = count_hydrogen(atoms)
                if "top" in h_folder:
                    gvib = gvib_source.get('top', 0.0) * n_H / 8
                elif "hol" in h_folder:
                    gvib = gvib_source.get('hol', 0.0) * n_H / 8
                else:
                    gvib = 0.0
                G0 = calculate_free_energy(energy, energy0, n_H, 0, gvib)
                data.append({
                    'label': get_label_from_subfolder(subfolder),
                    'energy': energy,
                    'energy0': energy0,
                    'n_H': n_H,
                    'G_vib': gvib,
                    'G0': G0
                })

            save_data_to_files(data, U_range, base_path, folder_stripped, h_folder_stripped)

            h_ads_results = calc_h_adsorption_energies(data, g_h2)
            print(f"\n[{folder_stripped}/{h_folder_stripped}] H adsorption ΔG (U=0):")
            for res in h_ads_results:
                print(f"{res['from']} → {res['to']}: ΔG_H = {res['delta_G_H']:.4f} eV")

            plt.figure(figsize=(5, 4))
            plt.plot(U_range, [0, 0], color='black', label='clean')
            for entry in data:
                Gmin = calculate_free_energy(entry['energy'], energy0, entry['n_H'], Umin, entry['G_vib'])
                Gmax = calculate_free_energy(entry['energy'], energy0, entry['n_H'], Umax, entry['G_vib'])
                plt.plot(U_range, [Gmin, Gmax], label=entry['label'])

            plt.xlabel('Potential (V vs RHE)')
            plt.ylabel('Relative Gibbs Free Energy (eV)')
            plt.xlim(Umin, Umax)
            plt.ylim(ymin, ymax)
            plt.legend(bbox_to_anchor=(1.05, 1.05), loc='upper left', fontsize='small', labelspacing=0.2)
            save_path = os.path.join(base_path, f"hbe_{folder_stripped}_{h_folder_stripped}.png")
            plt.savefig(save_path, bbox_inches='tight')
            plt.close()
            print(f"Saved: {save_path}")

if __name__ == "__main__":
    main() 