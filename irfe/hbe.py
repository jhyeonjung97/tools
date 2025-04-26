import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
import pandas as pd
import re

# 전압 범위 정의
Umin, Umax = -0.5, 1.5
ymin, ymax = -2.0, 1.0

def count_elements(atoms):
    symbols = atoms.get_chemical_symbols()
    return {'O': symbols.count('O'), 'H': symbols.count('H')}

def calculate_free_energy(E, E0, n_H, n_O, U, G_vib=0.0, pH=0):
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
        G_U = [calculate_free_energy(entry['energy'], entry['E_clean'], entry['n_H'], entry['n_O'], U, entry['G_vib']) for U in U_range]
        row_data = {
            'label': entry['label'],
            'n_H': entry['n_H'],
            'n_O': entry['n_O'],
            'E': entry['energy'],
            'E0': entry['E_clean'],
            'G_vib': entry['G_vib'],
            f'G(U={U_range[0]:.2f})': G_U[0],
            f'G(U={U_range[1]:.2f})': G_U[1]
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

def main():
    base_path = "/Users/hailey/Desktop/4_IrFe3"
    U_range = [Umin, Umax]

    # G_vib 값 미리 읽기
    gvib_ir = get_gvib_dict(os.path.join(base_path, "pourbaix_data_Ir.csv"))
    gvib_irfe = get_gvib_dict(os.path.join(base_path, "pourbaix_data_IrFe.csv"))

    for folder in sorted(os.listdir(base_path)):
        if not ("Ir" in folder and os.path.isdir(os.path.join(base_path, folder))):
            continue
        folder_path = os.path.join(base_path, folder)
        folder_stripped = strip_leading_number(folder)
        # *_H_* 폴더만 선택
        for h_folder in sorted(os.listdir(folder_path)):
            if not ("_H_" in h_folder and os.path.isdir(os.path.join(folder_path, h_folder))):
                continue
            h_folder_path = os.path.join(folder_path, h_folder)
            h_folder_stripped = strip_leading_number(h_folder)
            # clean surface: 1H의 final_with_calculator.json을 기준으로 사용
            clean_path = os.path.join(h_folder_path, "1H", "final_with_calculator.json")
            if not os.path.exists(clean_path):
                print(f"Clean surface not found: {clean_path}")
                continue
            clean_atoms = read(clean_path)
            E_clean = clean_atoms.get_potential_energy()

            # 어떤 시스템인지에 따라 G_vib 소스 선택
            if "IrFe" in folder:
                gvib_source = gvib_irfe
            else:
                gvib_source = gvib_ir

            data = []
            for subfolder in sorted(os.listdir(h_folder_path)):
                subfolder_path = os.path.join(h_folder_path, subfolder)
                json_path = os.path.join(subfolder_path, "final_with_calculator.json")
                if not os.path.exists(json_path):
                    continue
                atoms = read(json_path)
                energy = atoms.get_potential_energy()
                counts = count_elements(atoms)
                n_H = counts['H']
                if "top" in h_folder:
                    gvib = gvib_source.get('top', 0.0) * n_H / 8
                elif "hol" in h_folder:
                    gvib = gvib_source.get('hol', 0.0) * n_H / 8
                else:
                    gvib = 0.0
                data.append({
                    'label': get_label_from_subfolder(subfolder),
                    'energy': energy,
                    'E_clean': E_clean,
                    'n_H': n_H,
                    'n_O': counts['O'],
                    'G_vib': gvib
                })

            # 데이터 저장
            save_data_to_files(data, U_range, base_path, folder_stripped, h_folder_stripped)

            # 그래프 그리기
            plt.figure(figsize=(5, 4))
            plt.plot(U_range, [0, 0], color='black', label='clean')
            for entry in data:
                Gmin = calculate_free_energy(entry['energy'], E_clean, entry['n_H'], entry['n_O'], Umin, entry['G_vib'])
                Gmax = calculate_free_energy(entry['energy'], E_clean, entry['n_H'], entry['n_O'], Umax, entry['G_vib'])
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