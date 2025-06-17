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
    selected_systems = ['1_Mn', '2_Fe', '3_Co', '4_Ni', '0_Ir']
    system_names = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']
    
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
        json_path = os.path.join(system_dir, "final.json")
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
        if system == '0_Ir':
            selected_sites = ['1_layer_top', '2_layer_hol', '3_vac_top', '4_vac_hol', '5_atom_top', '6_atom_hol']
            site_names = ['ll_top', 'll_hol', 'vac_top', 'vac_hol', 'Ir_top', 'Ir_hol']
        else:
            selected_sites = ['2_layer_hol', '3_vac_hol', '4_vac_hol', '5_Ir_top', '6_Ir_hol', '8_M_hol']
            site_names = ['ll_hol', 'vac_hol', 'vac_hol', 'Ir_top', 'Ir_hol', 'M_hol']
                
        for target_site in selected_sites:

            # 지정된 사이트만 처리
            site_path = os.path.join(system_dir, target_site)
            if not os.path.exists(site_path) or not os.path.isdir(site_path):
                print(f"Warning: Target site {target_site} not found for {system}")
                continue
            
            # DONE 파일 확인
            if not is_calculation_done(site_path):
                print(f"Warning: DONE file not found for {system}/{target_site}, skipping...")
                continue

            # vib.txt 파일 확인
            vib_path = os.path.join(site_path, "vib.txt")
            if not os.path.exists(vib_path):
                print(f"Warning: No vib.txt for {system}/{target_site}")
                continue

            # final_with_calculator.json 파일 확인
            json_path = os.path.join(site_path, "final.json")
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
                if 'layer' in target_site.lower():
                    n_hydrogen = 8
                elif 'vac' in target_site.lower():
                    n_hydrogen = 7
                else:
                    n_hydrogen = 1
                
                # 흡착 에너지 계산
                binding_energy = calculate_adsorption_energy(
                    ads_energy, 
                    clean_energy, 
                    n_hydrogen,
                    h2_energy,
                    vib_correction
                )
                
                # 결과 저장
                system_name = system_names[selected_systems.index(system)]
                site_name = site_names[selected_sites.index(target_site)]
                result = {
                    'System': system_name,
                    'Site': site_name,
                    'Clean_energy': clean_energy,
                    'Ads_energy': ads_energy,
                    'H_count': n_hydrogen,
                    'G_vib': vib_correction,
                    'Binding_energy': binding_energy
                }
                
                results.append(result)
                print(f"  {system}/{target_site}: {n_hydrogen} H atoms, {binding_energy:.4f} eV, vib: {vib_correction:.4f} eV")
            except Exception as e:
                print(f"Error processing {json_path}: {e}")
    
    for surf in ['IrMn', 'IrFe', 'IrCo', 'IrNi']:
        for result in results:
            if result['System'] == surf and 'vac' in result['Site']:
                hbe_8 = next(r['Binding_energy'] for r in results if r['System'] == surf and r['H_count'] == 8)
                hbe_7 = result['Binding_energy']
                result['Binding_energy'] = hbe_8*8 - hbe_7*7

    for site in ['top', 'hol']:
        hbe_8 = next(r['Binding_energy'] for r in results if r['System'] == 'Ir' and r['Site'] == f'll_{site}')
        hbe_7 = next(r['Binding_energy'] for r in results if r['System'] == 'Ir' and r['Site'] == f'vac_{site}')
        for result in results:
            if result['System'] == 'Ir' and result['Site'] == f'vac_{site}':
                result['Binding_energy'] = hbe_8*8 - hbe_7*7

    # 결과를 DataFrame으로 변환    
    if results:
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
    수소 흡착 에너지 그래프 생성
    x축: 표면 종류
    y축: 흡착 에너지
    각 site마다 다른 색상과 마커로 표시
    """
    plt.figure(figsize=(10, 6))
    
    # 시스템 순서와 이름 정의
    system_order = {'IrMn': 0, 'IrFe': 1, 'IrCo': 2, 'IrNi': 3, 'Ir': 4}
    df['Order'] = df['System'].map(system_order)
    df = df.sort_values('Order')
    
    # 각 site별로 다른 색상과 마커 설정
    sites = df['Site'].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(sites)))  # Set2 컬러맵 사용
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']  # 다양한 마커
    
    # 각 site별로 데이터 플롯
    for i, site in enumerate(sites):
        site_data = df[df['Site'] == site]
        plt.scatter(site_data['Order'], 
                   site_data['Binding_energy'],
                   color=colors[i],
                   marker=markers[i % len(markers)],
                   label=site,
                   s=100,  # 마커 크기
                   alpha=0.7)
    
    # 이상적인 흡착 에너지 값 (약 0 eV)
    plt.axhline(y=0, color='black', linestyle='--', zorder=1, linewidth=1.0)
    
    # x축 레이블 설정
    x_ticks = range(len(df['System'].unique()))
    x_labels = sorted(df['System'].unique(), key=lambda x: system_order[x])
    plt.xticks(x_ticks, x_labels, rotation=45)
    
    # 그래프 꾸미기
    plt.ylabel('ΔG$_{\mathrm{H}}$ (eV)')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 범례 추가
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, framealpha=1.0)
    
    # x축 범위 설정 (약간의 여백)
    plt.xlim(-0.5, len(x_labels) - 0.5)
    
    # 레이아웃 조정
    plt.tight_layout()
    
    # 저장 및 표시
    plot_output = os.path.join(base_dir, 'H_binding_energies_final.png')
    plt.savefig(plot_output, bbox_inches='tight', dpi=300)
    print(f"Plot saved to {plot_output}")
    plt.close()

if __name__ == "__main__":
    main() 