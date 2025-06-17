import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
import pandas as pd

def count_adsorbates(atoms, adsorbate_type='H'):
    """원자 객체에서 흡착물 원자 개수 계산"""
    symbols = atoms.get_chemical_symbols()
    return symbols.count('H')

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

def calculate_adsorption_energy(ads_energy, clean_energy, n_adsorbate=1, ref_energy=None, vib_correction=0.0):
    """
    H 흡착 에너지 계산 함수
    - H: E_bind = (E_ads - E_clean - n_H * 0.5 * E_H2) / n_H
    """
    if ref_energy is None:
        ref_energy = -6.77108058  # H2 분자의 DFT 에너지
    return (ads_energy - clean_energy - n_adsorbate * 0.5 * ref_energy + vib_correction) / n_adsorbate

def is_calculation_done(directory):
    """DONE 파일이 디렉토리에 존재하는지 확인"""
    done_file = os.path.join(directory, "DONE")
    return os.path.exists(done_file)

def main():
    # 기본 디렉토리 설정
    base_dir = "/Users/hailey/Desktop/4_IrFe3"
    
    # H2 분자 에너지 설정
    h2_energy = -6.77108058
    
    # 사용할 표면 리스트
    selected_systems = ['5_IrMn', '6_IrFe', '7_IrCo', '8_IrNi', '0_Ir']
    system_names = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']
    
    # 각 시스템별로 사용할 사이트 정의
    target_sites = {
        '5_IrMn': '4_atom_top2',
        '6_IrFe': '4_atom_top2',
        '7_IrCo': '4_atom_top2',
        '8_IrNi': '4_atom_top2',
        '0_Ir': ['1_layer_top', '3_atom_top']  # Ir에 대해 두 개의 사이트 사용
    }
    
    # 결과 저장을 위한 리스트
    results = []
    
    # 1. Clean surface 에너지 읽기
    clean_energies = {}
    slab_dir = os.path.join(base_dir, "slab")
    
    print("Reading clean surface energies...")
    for system in selected_systems:
        system_dir = os.path.join(slab_dir, system)
        
        # DONE 파일 확인
        if not is_calculation_done(system_dir):
            print(f"Warning: DONE file not found for {system}, skipping...")
            continue
        
        # final_with_calculator.json 파일 확인
        json_path = os.path.join(system_dir, "final_with_calculator.json")
        if not os.path.exists(json_path):
            print(f"Warning: Clean surface file not found for {system}")
            continue
        
        # ASE로 파일 읽기
        try:
            atoms = read(json_path)
            energy = atoms.get_potential_energy()
            clean_energies[system] = energy
            print(f"  {system}: {energy:.4f} eV")
        except Exception as e:
            print(f"Error reading {json_path}: {e}")
    
    # 2. 수소 흡착 에너지 계산 (1_H 디렉토리만 처리)
    ads_dir = '1_H'
    ads_path = os.path.join(base_dir, ads_dir)
    print(f"\nReading H adsorption data...")
    
    for system in selected_systems:
        system_dir = os.path.join(ads_path, system)
        if not os.path.exists(system_dir) or not os.path.isdir(system_dir):
            print(f"Warning: H adsorption directory not found for {system}")
            continue
        
        # Clean surface 에너지 확인
        if system not in clean_energies:
            print(f"Warning: No clean surface energy for {system}, skipping...")
            continue
        
        clean_energy = clean_energies[system]
        target_site_list = target_sites[system] if isinstance(target_sites[system], list) else [target_sites[system]]
        
        for target_site in target_site_list:
            # 지정된 사이트만 처리
            site_path = os.path.join(system_dir, target_site)
            if not os.path.exists(site_path) or not os.path.isdir(site_path):
                print(f"Warning: Target site {target_site} not found for {system}")
                continue
            
            # DONE 파일 확인
            if not is_calculation_done(site_path):
                print(f"Warning: DONE file not found for {system}/{target_site}, skipping...")
                continue
            
            # final_with_calculator.json 파일 확인
            json_path = os.path.join(site_path, "final_with_calculator.json")
            if not os.path.exists(json_path):
                print(f"Warning: No final data for {system}/{target_site}")
                continue
            
            # 진동 자유에너지 보정값 읽기
            vib_file = os.path.join(site_path, "vib.txt")
            vib_correction = read_vib_correction(vib_file)
            
            # ASE로 파일 읽기
            try:
                atoms = read(json_path)
                ads_energy = atoms.get_potential_energy()
                
                # 수소 원자 개수 계산
                n_hydrogen = count_adsorbates(atoms)
                
                # 레이어 여부에 따른 normalization 처리
                is_layer = 'layer' in target_site.lower()
                
                # 흡착 에너지 계산
                binding_energy = calculate_adsorption_energy(
                    ads_energy, clean_energy, 
                    n_hydrogen if is_layer else 1, 
                    h2_energy,
                    vib_correction
                )
                
                # 결과 저장
                system_name = system_names[selected_systems.index(system)]
                result = {
                    'System': system,
                    'System_name': system_name,
                    'Site': target_site,
                    'Clean_energy': clean_energy,
                    'Ads_energy': ads_energy,
                    'H_count': n_hydrogen,
                    'Vib_correction': vib_correction,
                    'Binding_energy': binding_energy
                }
                
                results.append(result)
                print(f"  {system}/{target_site}: {n_hydrogen} H atoms, {binding_energy:.4f} eV{'(normalized)' if is_layer else ''}, vib: {vib_correction:.4f} eV")
            except Exception as e:
                print(f"Error processing {json_path}: {e}")
    
    # 결과를 DataFrame으로 변환
    if results:
        print(results)
        df = pd.DataFrame(results)
        
        # CSV 및 TSV로 저장
        csv_output = os.path.join(base_dir, 'H_binding_energies_final.csv')
        tsv_output = os.path.join(base_dir, 'H_binding_energies_final.tsv')
        
        df.to_csv(csv_output, index=False)
        df.to_csv(tsv_output, sep='\t', index=False, float_format='%.4f')
        
        print(f"Results saved to {csv_output} and {tsv_output}")
        
        # 그래프 그리기
        plot_binding_energies(df, base_dir)

def plot_binding_energies(df, base_dir):
    """
    수소 흡착 에너지 그래프 생성 (간략화된 버전)
    x축: 표면 종류
    y축: 흡착 에너지
    """
    plt.figure(figsize=(4, 3))
    
    # 시스템 순서와 이름 정의
    system_order = {'5_IrMn': 0, '6_IrFe': 1, '7_IrCo': 2, '8_IrNi': 3, '0_Ir': 4}
    df['Order'] = df['System'].map(system_order)
    df = df.sort_values('Order')
    
    # 데이터 플롯
    x = df['Order']
    y = df['Binding_energy']
    
    # 4_atom_top2 데이터 포인트 플롯 (Ir의 3_atom_top 포함)
    top2_mask = ((df['System'] != '0_Ir') & (df['Site'] == '4_atom_top2')) | \
                ((df['System'] == '0_Ir') & (df['Site'] == '3_atom_top'))
    plt.plot(x[top2_mask], y[top2_mask], color='black', marker='o', linestyle='-', 
             label='low H coverage', zorder=2)
    
    # Ir의 1_layer_top 데이터 포인트 플롯
    ir_layer_mask = (df['System'] == '0_Ir') & (df['Site'] == '1_layer_top')
    plt.plot(x[ir_layer_mask], y[ir_layer_mask], color='black', marker='o', markerfacecolor='white', 
             markeredgecolor='black', linestyle='-', label='high H coverage', zorder=2)
    
    # 이상적인 흡착 에너지 값 (약 0 eV)
    plt.axhline(y=0, color='black', linestyle='--', zorder=1, linewidth=1.0)
    
    # x축 레이블 설정
    x_ticks = range(len(df['System_name'].unique()))
    x_labels = [name for name in df.sort_values('Order')['System_name'].unique()]
    plt.xticks(x_ticks, x_labels)
    
    # 그래프 꾸미기
    plt.ylabel('ΔG$_{\mathrm{H}}$ (eV)')
    
    # 범례 추가
    plt.legend(loc='upper right', frameon=True, framealpha=1.0)
    
    # x축 범위 설정 (약간의 여백)
    plt.xlim(-0.5, len(x_labels) - 0.5)
    # plt.ylim(-0.45, 0.25)
    
    # 저장 및 표시
    plot_output = os.path.join(base_dir, 'H_binding_energies_final.png')
    plt.savefig(plot_output, bbox_inches='tight')
    print(f"Plot saved to {plot_output}")
    plt.close()

if __name__ == "__main__":
    main() 