import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
import pandas as pd
import re

def count_adsorbates(atoms, adsorbate_type='H'):
    """원자 객체에서 흡착물 원자 개수 계산"""
    symbols = atoms.get_chemical_symbols()
    if adsorbate_type == 'H':
        return symbols.count('H')
    elif adsorbate_type == 'OH':
        # OH는 O와 H의 쌍으로 구성
        return symbols.count('O')  # O 개수로 OH 개수 계산
    elif adsorbate_type == 'O':
        return symbols.count('O')
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

def calculate_adsorption_energy(ads_energy, clean_energy, n_adsorbate=1, ref_energy=None, adsorbate_type='H', vib_correction=0.0):
    """
    흡착 에너지 계산 함수
    - H: E_bind = (E_ads - E_clean - n_H * 0.5 * E_H2) / n_H
    - OH: E_bind = (E_ads - E_clean - n_OH * (E_H2O - 0.5 * E_H2)) / n_OH
    - O: E_bind = (E_ads - E_clean - n_O * 0.5 * E_O2) / n_O
    """
    if adsorbate_type == 'H':
        # H2 분자의 DFT 에너지
        if ref_energy is None:
            ref_energy = -6.77108058
        return (ads_energy - clean_energy - n_adsorbate * 0.5 * ref_energy + vib_correction) / n_adsorbate
    
    elif adsorbate_type == 'OH':
        # H2O 및 H2 분자의 DFT 에너지
        if ref_energy is None:
            h2o_energy = -14.22930  # H2O 에너지
            h2_energy = -6.77108058  # H2 에너지
            ref_energy = h2o_energy - 0.5 * h2_energy  # OH 레퍼런스
        return (ads_energy - clean_energy - n_adsorbate * ref_energy + vib_correction) / n_adsorbate
    
    elif adsorbate_type == 'O':
        # O2 분자의 DFT 에너지
        if ref_energy is None:
            ref_energy = -9.85663  # O2 에너지
        return (ads_energy - clean_energy - n_adsorbate * 0.5 * ref_energy + vib_correction) / n_adsorbate
    
    return 0

def is_calculation_done(directory):
    """DONE 파일이 디렉토리에 존재하는지 확인"""
    done_file = os.path.join(directory, "DONE")
    return os.path.exists(done_file)

def main():
    # 기본 디렉토리 설정
    base_dir = "/Users/hailey/Desktop/4_IrFe3"
    
    # 레퍼런스 에너지 설정
    h2_energy = -6.77108058   # H2 분자 에너지
    o2_energy = -9.85663      # O2 분자 에너지
    h2o_energy = -14.22930    # H2O 분자 에너지
    oh_ref_energy = h2o_energy - 0.5 * h2_energy  # OH 레퍼런스 에너지
    
    # 사용할 표면 리스트
    selected_systems = ['5_IrMn', '6_IrFe', '7_IrCo', '8_IrNi', '0_Ir']
    system_names = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']
    
    # 흡착물 리스트
    adsorbates = {
        '1_H': {'name': 'H', 'ref_energy': h2_energy, 'type': 'H'},
        '2_OH': {'name': 'OH', 'ref_energy': oh_ref_energy, 'type': 'OH'},
        '3_O': {'name': 'O', 'ref_energy': o2_energy, 'type': 'O'}
    }
    
    # 결과 저장을 위한 리스트
    all_results = []
    
    # 1. Clean surface 에너지 읽기
    clean_energies = {}
    clean_atoms = {}  # Clean surface 원자 객체도 저장
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
            clean_atoms[system] = atoms
            print(f"  {system}: {energy:.4f} eV")
        except Exception as e:
            print(f"Error reading {json_path}: {e}")
    
    # 2. 각 흡착물에 대한 흡착 에너지 계산
    for ads_dir, ads_info in adsorbates.items():
        ads_name = ads_info['name']
        ads_ref_energy = ads_info['ref_energy']
        ads_type = ads_info['type']
        
        # 결과 저장을 위한 딕셔너리
        results = []
        
        ads_path = os.path.join(base_dir, ads_dir)
        print(f"\nReading {ads_name} adsorption data...")
        
        for system in selected_systems:
            system_dir = os.path.join(ads_path, system)
            if not os.path.exists(system_dir) or not os.path.isdir(system_dir):
                print(f"Warning: {ads_name} adsorption directory not found for {system}")
                continue
            
            # Clean surface 에너지 확인
            if system not in clean_energies:
                print(f"Warning: No clean surface energy for {system}, skipping...")
                continue
            
            clean_energy = clean_energies[system]
            
            # Coverage_site 서브폴더 확인
            for site_dir in os.listdir(system_dir):
                site_path = os.path.join(system_dir, site_dir)
                if not os.path.isdir(site_path):
                    continue
                
                # DONE 파일 확인
                if not is_calculation_done(site_path):
                    print(f"Warning: DONE file not found for {system}/{site_dir}, skipping...")
                    continue
                
                # final_with_calculator.json 파일 확인
                json_path = os.path.join(site_path, "final_with_calculator.json")
                if not os.path.exists(json_path):
                    print(f"Warning: No final data for {system}/{site_dir}")
                    continue
                
                # 진동 자유에너지 보정값 읽기
                vib_file = os.path.join(site_path, "vib", "vib.txt")
                vib_correction = read_vib_correction(vib_file)
                
                # ASE로 파일 읽기
                try:
                    atoms = read(json_path)
                    ads_energy = atoms.get_potential_energy()
                    
                    # 흡착물 개수 계산
                    n_adsorbate = count_adsorbates(atoms, ads_type)
                    
                    # Layer 여부에 따른 normalization 처리
                    is_layer = 'layer' in site_dir.lower()
                    
                    # 흡착 에너지 계산
                    binding_energy = calculate_adsorption_energy(
                        ads_energy, clean_energy, 
                        n_adsorbate if is_layer else 1, 
                        ads_ref_energy, ads_type,
                        vib_correction
                    )
                    
                    # 결과 저장
                    system_name = system_names[selected_systems.index(system)]
                    result = {
                        'Adsorbate': ads_name,
                        'System': system,
                        'System_name': system_name,
                        'Site': site_dir,
                        'Clean_energy': clean_energy,
                        'Ads_energy': ads_energy,
                        'Ads_count': n_adsorbate,
                        'Is_layer': is_layer,
                        'Vib_correction': vib_correction,
                        'Binding_energy': binding_energy
                    }
                    
                    results.append(result)
                    all_results.append(result)
                    
                    print(f"  {system}/{site_dir}: {n_adsorbate} {ads_name} atoms, {binding_energy:.4f} eV{'(normalized)' if is_layer else ''}, vib: {vib_correction:.4f} eV")
                except Exception as e:
                    print(f"Error processing {json_path}: {e}")
        
        # 결과를 DataFrame으로 변환
        if results:
            df = pd.DataFrame(results)
            
            # CSV 및 TSV로 저장
            csv_output = os.path.join(base_dir, 'figures', f'{ads_name}_binding_energies.csv')
            tsv_output = os.path.join(base_dir, 'figures', f'{ads_name}_binding_energies.tsv')
            
            df.to_csv(csv_output, index=False)
            df.to_csv(tsv_output, sep='\t', index=False, float_format='%.4f')
            
            print(f"Results saved to {csv_output} and {tsv_output}")
            
            # 그래프 그리기
            plot_binding_energies(df, base_dir, ads_name)
    
    # 모든 흡착물 데이터 저장
    if all_results:
        df_all = pd.DataFrame(all_results)
        csv_output = os.path.join(base_dir, 'figures', 'all_binding_energies.csv')
        tsv_output = os.path.join(base_dir, 'figures', 'all_binding_energies.tsv')
        
        df_all.to_csv(csv_output, index=False)
        df_all.to_csv(tsv_output, sep='\t', index=False, float_format='%.4f')
        
        print(f"\nAll results saved to {csv_output} and {tsv_output}")

def plot_binding_energies(df, base_dir, adsorbate='H'):
    """
    흡착 에너지 그래프 생성
    x축: 표면 종류
    y축: 흡착 에너지
    범례: coverage_site
    """
    plt.figure(figsize=(4, 3))
    
    # 시스템 순서와 이름 정의
    system_order = {'5_IrMn': 0, '6_IrFe': 1, '7_IrCo': 2, '8_IrNi': 3, '0_Ir': 4}
    df['Order'] = df['System'].map(system_order)
    df = df.sort_values('Order')
    
    # 고유한 사이트 목록 가져오기 (이름순 정렬)
    sites = sorted(df['Site'].unique())
    
    # 컬러맵 설정
    colors = plt.cm.tab10(np.linspace(0, 1, len(sites)))
    
    # 각 사이트별로 그래프 그리기
    for i, site in enumerate(sites):
        site_data = df[df['Site'] == site]
        x = site_data['Order']
        y = site_data['Binding_energy']
        
        plt.scatter(x, y, color=colors[i], label=f"{site} ({site_data['Ads_count'].iloc[0]} {adsorbate})", zorder=3)
        plt.plot(x, y, color=colors[i], linestyle='-', zorder=2)
    
    # 이상적인 흡착 에너지 값 (약 0 eV)
    plt.axhline(y=0, color='gray', linestyle='--', zorder=1)
    
    # x축 레이블 설정
    x_ticks = range(len(df['System_name'].unique()))
    x_labels = [name for name in df.sort_values('Order')['System_name'].unique()]
    plt.xticks(x_ticks, x_labels)
    
    # 그래프 꾸미기
    plt.ylabel(f'{adsorbate} Binding Energy (eV)')
    plt.legend(title='Coverage & Site', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # x축 범위 설정 (약간의 여백)
    plt.xlim(-0.5, len(x_labels) - 0.5)
    
    # 저장 및 표시
    plot_output = os.path.join(base_dir, 'figures', f'{adsorbate}_binding_energies.png')
    os.makedirs(os.path.dirname(plot_output), exist_ok=True)
    plt.savefig(plot_output, bbox_inches='tight')
    print(f"Plot saved to {plot_output}")
    plt.close()

if __name__ == "__main__":
    main() 