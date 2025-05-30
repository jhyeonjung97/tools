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

def get_site_number(folder_name):
    """사이트 위치에 따른 숫자 반환 (colormap 인덱스용)"""
    if 'layer_top' in folder_name:
        return 0
    elif 'layer_hol' in folder_name:
        return 1
    elif 'Ir_top' in folder_name:
        return 2
    elif 'Ir_hol' in folder_name:
        return 3
    elif 'M_top' in folder_name:
        return 4
    elif 'M_hol' in folder_name:
        return 5
    return 0

def read_vib_correction(vib_file):
    """vib.txt 파일에서 G(T) 보정값을 읽는 함수"""
    try:
        with open(vib_file, 'r') as f:
            content = f.read()
            # Thermal correction to G(T) 값을 찾아서 eV 단위로 반환
            for line in content.split('\n'):
                if 'Zero-point energy E_ZPE' in line:
                    zpe = float(line.split()[-2])  # eV 값 반환
                if 'Entropy contribution T*S' in line:
                    ts = float(line.split()[-2])  # eV 값 반환
            return zpe - ts
    except:
        return 0.0  # 파일이 없거나 읽을 수 없는 경우 0 반환

def calculate_free_energy(E, E0, n_H, n_O, U, G_vib=0.0, pH=0):
    """Gibbs free energy 계산 - 진동 자유에너지와 열역학적 보정 포함"""
    # 기체 상수
    kT = 0.026  # kT at room temperature (eV)
    const = kT * np.log(10)  # 0.0592 eV
    
    # 기체 분자의 에너지와 보정값
    e_h2 = -6.77108058  # H2의 free energy (eV)
    e_h2o = -14.22266334  # H2O의 free energy (eV)
    
    # ZPE, 열용량, 엔트로피 보정
    zpe_h2o = 0.558
    cv_h2o = 0.103
    ts_h2o = 0.675

    zpe_h2 = 0.268
    cv_h2 = 0.0905
    ts_h2 = 0.408

    # 기체 분자의 깁스 자유에너지
    g_h2o = e_h2o + zpe_h2o - ts_h2o + cv_h2o
    g_h2 = e_h2 + zpe_h2 - ts_h2 + cv_h2

    # 흡착물의 깁스 자유에너지 계산
    surface_term = E - E0 + G_vib  # 전자 에너지 + 진동 자유에너지
    
    # H*, OH*, O* 흡착물의 전하 상태에 따른 U 의존성
    # H*: +1e, OH*: -1e, O*: -2e, OOH*: -3e 기준
    n_e = (n_H - 2*n_O)  # H*는 +1e, O*와 OH*는 -e 기여
    
    # 전체 깁스 자유에너지 계산
    G = surface_term + n_e * U + n_e * const * pH
    
    # 기체상 보정
    G -= n_H * 0.5 * g_h2  # H2 contribution
    G -= n_O * (g_h2o - g_h2)  # H2O contribution
    
    return G/8

def clean_label(ads, site):
    """폴더명에서 숫자를 제거하고 깔끔한 레이블 생성"""
    # 흡착물 이름에서 숫자 제거
    clean_ads = ads.split('_')[1]
    # 사이트 이름에서 숫자 제거
    clean_site = '_'.join(s for s in site.split('_') if not s[0].isdigit())
    return f"{clean_ads}_{clean_site}"

def get_sort_key(entry):
    """정렬을 위한 키 생성"""
    # 흡착물 순서: H, OH, O
    ads_order = {'H': 0, 'OH': 1, 'O': 2}
    ads = entry['ads'].split('_')[1]
    
    # 사이트 순서: layer_top, layer_hol, Ir_top, Ir_hol, M_top, M_hol
    site_order = {
        'layer_top': 0,
        'layer_hol': 1,
        'Ir_top': 2,
        'Ir_hol': 3,
        'M_top': 4,
        'M_hol': 5
    }
    
    site = entry['site']
    for key in site_order.keys():
        if key in site:
            return (ads_order[ads], site_order[key])
    return (ads_order[ads], 999)  # 기타 경우

def save_data_to_files(data, base_path, system_name):
    """데이터를 CSV와 TSV 파일로 저장"""
    # 데이터 프레임 생성을 위한 리스트
    df_data = []
    
    for entry in data:
        # site 이름에서 숫자 제거하고 분리
        site_parts = [s for s in entry['site'].split('_') if not s[0].isdigit()]
        site1 = site_parts[0] if len(site_parts) > 0 else ''
        site2 = site_parts[1] if len(site_parts) > 1 else ''
        
        # 각 전압에 대한 에너지 값을 행으로 추가
        row_data = {
            'ads': entry['ads'],
            'site1': site1,
            'site2': site2,
            'E': entry['E'],
            'E0': entry['E0'],
            'n_H': entry['n_H'],
            'n_O': entry['n_O'],
            'G_vib': entry['G_vib'],
            'Gmin': entry['Gmin'],
            'Gmax': entry['Gmax'],
            'G0': entry['G0']
        }
        df_data.append(row_data)
    
    # DataFrame 생성
    df = pd.DataFrame(df_data)
    
    # CSV 파일 저장 (모든 소수점 유지)
    csv_path = os.path.join(base_path, 'figures', f'pourbaix_data_{system_name}.csv')
    df.to_csv(csv_path, index=False)
    
    # TSV 파일 저장 (소수점 2자리까지)
    tsv_path = os.path.join(base_path, 'figures', f'pourbaix_data_{system_name}.tsv')
    df_rounded = df.round(2)
    df_rounded.to_csv(tsv_path, sep='\t', index=False)

def create_pourbaix_diagram(base_path):
    # slab 하위 폴더명에서 metal_systems 자동 생성
    slab_path = os.path.join(base_path, 'slab')
    metal_systems = {}
    if os.path.exists(slab_path):
        for folder in os.listdir(slab_path):
            if not folder.startswith('.'):
                # 예: 0_Ir, 2_Fe 등 폴더명에서 시스템명 추출
                if '_' in folder:
                    sys_name = folder.split('_', 1)[1]
                else:
                    sys_name = folder
                metal_systems[folder] = sys_name

    # 색상 및 alpha 설정 딕셔너리
    color_alpha_dict = {
        'clean': ('black', 1.0),
        'H_layer': ('b', 1.0),
        'H_atom': ('b', 0.5),
        'OH_layer': ('g', 1.0),
        'OH_atom': ('g', 0.5),
        'O_layer': ('r', 1.0),
        'O_atom': ('r', 0.5)
    }

    U_range = [Umin, Umax]

    cross_potentials = []  # (system_name, x_cross) 저장용 리스트
    for target_metal in metal_systems.keys():
        clean_path = os.path.join(base_path, 'slab', target_metal, 'final_with_calculator.json')
        if not os.path.exists(clean_path):
            continue
        clean_done_path = os.path.join(base_path, 'slab', target_metal, 'DONE')
        if not os.path.exists(clean_done_path):
            continue
        clean_atoms = read(clean_path)
        energy_clean = clean_atoms.get_potential_energy()

        all_data = []  # 모든 site 데이터 저장용
        stable_sites = {}  # group_key별 가장 안정한 것만 저장

        for ads in ['1_H', '2_OH', '3_O']:
            metal_path = os.path.join(base_path, ads, target_metal)
            if not os.path.exists(metal_path):
                continue
            for site in os.listdir(metal_path):
                if site.startswith('.') or site == 'vib':
                    continue
                if os.path.exists(os.path.join(metal_path, site, 'unmatched')):
                    continue
                full_path = os.path.join(metal_path, site, 'final_with_calculator.json')
                if not os.path.exists(full_path):
                    continue
                done_path = os.path.join(metal_path, site, 'DONE')
                if not os.path.exists(done_path):
                    continue
                vib_file = os.path.join(metal_path, site, 'vib.txt')
                G_vib = read_vib_correction(vib_file)
                atoms = read(full_path)
                energy = atoms.get_potential_energy()
                counts = count_elements(atoms)
                Gmin = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umin, G_vib)
                Gmax = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umax, G_vib)
                G0 = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], 0, G_vib)
                # 그룹핑 키: 흡착물 종류(H, OH, O) + 흡착 위치(layer, atom)
                if 'layer' in site:
                    group_key = f"{ads.split('_')[1]}_layer"
                elif 'atom' in site:
                    group_key = f"{ads.split('_')[1]}_atom"
                else:
                    group_key = f"{ads.split('_')[1]}_other"
                # 모든 데이터 저장
                all_data.append({
                    'ads': ads.split('_')[1],
                    'site': site,
                    'E': energy,
                    'E0': energy_clean,
                    'n_H': counts['H'],
                    'n_O': counts['O'],
                    'G_vib': G_vib,
                    'Gmin': Gmin,
                    'Gmax': Gmax,
                    'G0': G0
                })
                # 가장 안정한 것만 저장
                if group_key not in stable_sites or G0 < stable_sites[group_key]['G0']:
                    stable_sites[group_key] = {
                        'Gmin': Gmin,
                        'Gmax': Gmax,
                        'G0': G0,
                        'group_key': group_key
                    }

        # 모든 site 데이터 저장 (TSV/CSV)
        save_data_to_files(all_data, base_path, metal_systems[target_metal])

        # 그래프 생성 및 저장 (가장 안정한 것만)
        plt.figure(figsize=(5, 4))
        plt.plot(U_range, [0]*len(U_range), color='black', label='clean')
        # x=0 수직선 추가
        plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

        for group_key, entry in stable_sites.items():
            color, alpha = color_alpha_dict.get(group_key, ('gray', 1.0))
            plt.plot(U_range, [entry['Gmin'], entry['Gmax']], color=color, alpha=alpha, label=group_key)
        plt.xlabel('Potential (V vs RHE)')
        plt.ylabel('Relative Gibbs Free Energy (eV)')
        plt.xlim(Umin, Umax)
        plt.ylim(ymin, ymax)
        plt.legend(bbox_to_anchor=(1.05, 1.05), loc='upper left', fontsize='small', labelspacing=0.2)
        save_path = os.path.join(base_path, 'figures', f'pourbaix_diagram_{metal_systems[target_metal]}.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

        # OH_layer와 O_layer 교차점 포텐셜을 리스트에 저장
        if 'OH_layer' in stable_sites and 'O_layer' in stable_sites:
            oh = stable_sites['OH_layer']
            o = stable_sites['O_layer']
            a_oh = (oh['Gmax'] - oh['Gmin']) / (Umax - Umin)
            b_oh = oh['Gmin'] - a_oh * Umin
            a_o = (o['Gmax'] - o['Gmin']) / (Umax - Umin)
            b_o = o['Gmin'] - a_o * Umin
            if a_oh != a_o:
                x_cross = (b_o - b_oh) / (a_oh - a_o)
                if Umin <= x_cross <= Umax:
                    cross_potentials.append((metal_systems[target_metal], x_cross))

    # 원하는 순서대로 출력
    order = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']
    for sys in order:
        for name, x_cross in cross_potentials:
            if name == sys:
                print(f"{name}: OH_layer -> O_layer = {x_cross:.2f} V")

if __name__ == "__main__":
    base_path = "/Users/hailey/Desktop/4_IrFe3"
    create_pourbaix_diagram(base_path)
