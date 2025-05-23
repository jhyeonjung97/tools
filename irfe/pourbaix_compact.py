import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase import Atoms
import matplotlib.colors as mcolors
import pandas as pd

# 전압 범위 정의
Umin, Umax = -0.5, 2.0
ymin, ymax = -1.5, 0.5

def count_elements(atoms):
    """원자 구조에서 O와 H의 개수를 세는 함수"""
    symbols = atoms.get_chemical_symbols()
    return {'O': symbols.count('O'), 'H': symbols.count('H')}

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

def clean_label(ads):
    """흡착물 라벨 정리 (1_H -> H)"""
    if '_' in ads:
        return ads.split('_')[1]
    return ads

def save_data_to_files(data, base_path, system_name):
    """데이터를 CSV와 TSV 파일로 저장"""
    # 데이터 프레임 생성을 위한 리스트
    df_data = []
    
    for entry in data:
        # 각 전압에 대한 에너지 값을 행으로 추가
        row_data = {
            'ads': entry['ads'],
            'site': entry['site'],
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
    csv_path = os.path.join(base_path, f'pourbaix_data_{system_name}.csv')
    df.to_csv(csv_path, index=False)
    
    # TSV 파일 저장 (소수점 2자리까지)
    tsv_path = os.path.join(base_path, f'pourbaix_data_{system_name}.tsv')
    df_rounded = df.round(2)
    df_rounded.to_csv(tsv_path, sep='\t', index=False)

def is_calculation_done(directory):
    """DONE 파일이 디렉토리에 존재하는지 확인"""
    done_file = os.path.join(directory, "DONE")
    return os.path.exists(done_file)

def create_pourbaix_diagram(base_path):
    # 사용할 표면 리스트
    selected_systems = ['5_IrMn', '6_IrFe', '7_IrCo', '8_IrNi', '0_Ir']
    system_names = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']
    
    # 사용할 흡착물 리스트
    selected_adsorbates = ['1_H', '2_OH', '3_O']
    
    # 색상 설정 딕셔너리
    color_dict = {
        'H': 'blue',
        'OH': 'green',
        'O': 'red'
    }

    U_range = [Umin, Umax]

    cross_potentials = []  # (system_name, x_cross) 저장용 리스트
    
    for metal_code, system_name in zip(selected_systems, system_names):
        clean_path = os.path.join(base_path, 'slab', metal_code, 'final_with_calculator.json')
        if not os.path.exists(clean_path):
            continue
            
        # DONE 파일 확인
        if not is_calculation_done(os.path.join(base_path, 'slab', metal_code)):
            continue
            
        clean_atoms = read(clean_path)
        energy_clean = clean_atoms.get_potential_energy()

        all_data = []  # 모든 layer_top 데이터 저장용
        ads_data = {}  # 흡착물별 가장 안정한 layer_top 데이터

        for ads in selected_adsorbates:
            metal_path = os.path.join(base_path, ads, metal_code)
            if not os.path.exists(metal_path):
                continue
                
            for site in os.listdir(metal_path):
                # layer_top 데이터만 처리
                if 'layer_top' not in site:
                    continue
                    
                if site.startswith('.'):
                    continue
                    
                site_path = os.path.join(metal_path, site)
                if not os.path.isdir(site_path):
                    continue
                    
                # DONE 파일 확인
                if not is_calculation_done(site_path):
                    continue
                    
                full_path = os.path.join(site_path, 'final_with_calculator.json')
                if not os.path.exists(full_path):
                    continue
                    
                vib_file = os.path.join(site_path, 'vib.txt')
                G_vib = read_vib_correction(vib_file)
                atoms = read(full_path)
                energy = atoms.get_potential_energy()
                counts = count_elements(atoms)
                
                Gmin = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umin, G_vib)
                Gmax = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], Umax, G_vib)
                G0 = calculate_free_energy(energy, energy_clean, counts['H'], counts['O'], 0, G_vib)
                
                # 흡착물 이름 추출 (1_H -> H)
                ads_name = clean_label(ads)
                
                # 모든 데이터 저장
                all_data.append({
                    'ads': ads_name,
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
                if ads_name not in ads_data or G0 < ads_data[ads_name]['G0']:
                    ads_data[ads_name] = {
                        'Gmin': Gmin,
                        'Gmax': Gmax,
                        'G0': G0,
                        'site': site
                    }

        # 데이터가 없으면 건너뜀
        if not all_data:
            continue

        # 모든 데이터 저장 (TSV/CSV)
        save_data_to_files(all_data, base_path, system_name)

        # 그래프 생성 및 저장 (가장 안정한 것만)
        plt.figure(figsize=(4, 3))
        plt.plot(U_range, [0, 0], color='black', label='clean')
        # x=0 수직선 추가
        plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

        for ads_name, entry in ads_data.items():
            color = color_dict.get(ads_name, 'gray')
            plt.plot(U_range, [entry['Gmin'], entry['Gmax']], color=color, label=f"{ads_name}")
        
        plt.xlabel('Potential (V vs RHE)')
        plt.ylabel('ΔG (eV)')
        plt.xlim(Umin, Umax)
        plt.ylim(ymin, ymax)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plot_output = os.path.join(base_path, f'pourbaix_diagram_{system_name}.png')
        plt.savefig(plot_output, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to {plot_output}")

        # OH와 O 교차점 포텐셜을 리스트에 저장
        if 'OH' in ads_data and 'O' in ads_data:
            oh = ads_data['OH']
            o = ads_data['O']
            a_oh = (oh['Gmax'] - oh['Gmin']) / (Umax - Umin)
            b_oh = oh['Gmin'] - a_oh * Umin
            a_o = (o['Gmax'] - o['Gmin']) / (Umax - Umin)
            b_o = o['Gmin'] - a_o * Umin
            if a_oh != a_o:
                x_cross = (b_o - b_oh) / (a_oh - a_o)
                if Umin <= x_cross <= Umax:
                    cross_potentials.append((system_name, x_cross))

    # 원하는 순서대로 출력
    for sys in system_names:
        for name, x_cross in cross_potentials:
            if name == sys:
                print(f"{name}: OH -> O = {x_cross:.2f} V")

if __name__ == "__main__":
    base_path = "/Users/hailey/Desktop/4_IrFe3"
    create_pourbaix_diagram(base_path) 