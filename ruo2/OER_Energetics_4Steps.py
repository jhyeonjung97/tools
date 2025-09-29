from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# 사용자 환경에 맞춰 경로를 변경하세요
root = "~/Desktop/3_RuO2"

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

def get_energy_from_json(json_path):
    """JSON 파일에서 에너지를 읽어오는 함수"""
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def calculate_oer_energies():
    """14가지 OER 경로의 반응 에너지를 계산하고 표 형태로 정리하는 함수"""
    root_path = Path(root).expanduser()

    gibbs_correction_o_v = 0.058092 - 0.033706
    gibbs_correction_o_oh = 0.439910 - 0.128005
    gibbs_correction_o_o = 0.156862 - 0.100681
    gibbs_correction_oh_oo = 0.502349 - 0.144209

    # 표 형태로 데이터를 저장할 리스트
    oer_data = []
    all_paths = {}

    # Path 1: RuO2에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 1: RuO2 (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "2_RuO2_OER" / "7_RuO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path1'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 1: Failed to load some energies")
    except Exception as e:
        print(f"Path 1: Error - {e}")

    # Path 2: ReRuO2에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 2: ReRuO2 (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "3_ReRuO2_OER" / "3_ReRuO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path2'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 2: Failed to load some energies")
    except Exception as e:
        print(f"Path 2: Error - {e}")

    # Path 3: RuO2_strain_XRD에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 3: RuO2_strain_XRD (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "3_ReRuO2_OER" / "4_RuO2_strain_XRD"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2_1%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path3'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 3: Failed to load some energies")
    except Exception as e:
        print(f"Path 3: Error - {e}")

    # Path 4: RuO2_strain_ReO2에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 4: RuO2_strain_ReO2 (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "3_ReRuO2_OER" / "5_RuO2_strain_ReO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2_2%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path4'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 4: Failed to load some energies")
    except Exception as e:
        print(f"Path 4: Error - {e}")

    # Path 5: ReRuO2_strain_1%에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 5: ReRuO2_strain_1% (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "3_ReRuO2_OER" / "6_ReRuO2_strain_1%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2_1%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path5'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 5: Failed to load some energies")
    except Exception as e:
        print(f"Path 5: Error - {e}")

    # Path 6: ReRuO2_strain_2%에서 1_O_V → 2_O_OH → 3_O_O → 5_OH_OO
    print("\n=== Path 6: ReRuO2_strain_2% (O_V → O_OH → O_O → OH_OO) ===")
    try:
        oer_path = root_path / "3_ReRuO2_OER" / "7_ReRuO2_strain_2%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.json") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.json") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.json") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "4_OH_OO" / "final_with_calculator.json") + gibbs_correction_oh_oo
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step4 = 4.92 - step1 - step2 - step3
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2_2%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4
            })
            all_paths['Path6'] = [step1, step2, step3, step4]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV")
        else:
            print("Path 6: Failed to load some energies")
    except Exception as e:
        print(f"Path 6: Error - {e}")

    return all_paths, oer_data

def print_oer_table(oer_data):
    """OER 데이터를 표 형태로 출력하는 함수"""
    if not oer_data:
        print("No OER data to display")
        return
    
    print("\n" + "="*120)
    print("OER Energetics Summary Table")
    print("="*120)
    
    # 헤더 출력
    header = f"{'Surface':<12} {'Int1':<8} {'Int2':<8} {'Int3':<8} {'Int4':<8} {'Step1':<10} {'Step2':<10} {'Step3':<10} {'Step4':<10} {'Max':<10} {'Overpotential':<12}"
    print(header)
    print("-"*120)
    
    # 데이터 출력
    for i, data in enumerate(oer_data, 1):
        max_energy = max(data['step1'], data['step2'], data['step3'], data['step4'])
        overpotential = max_energy - 1.23
        
        row = f"{data['surface']:<12} {data['int1']:<8} {data['int2']:<8} {data['int3']:<8} {data['int4']:<8} " \
              f"{data['step1']:<10.3f} {data['step2']:<10.3f} {data['step3']:<10.3f} {data['step4']:<10.3f} " \
              f"{max_energy:<10.3f} {overpotential:<12.3f}"
        print(row)
    
    print("-"*120)
    
    # 최적 경로 찾기
    best_path = min(oer_data, key=lambda x: max(x['step1'], x['step2'], x['step3'], x['step4']))
    best_max = max(best_path['step1'], best_path['step2'], best_path['step3'], best_path['step4'])
    best_overpotential = best_max - 1.23
    
    print(f"\nBest Path: {best_path['surface']} - {best_path['int1']} → {best_path['int2']} → {best_path['int3']} → {best_path['int4']}")
    print(f"Maximum Step Energy: {best_max:.3f} eV")
    print(f"Overpotential: {best_overpotential:.3f} eV")
    print("="*120)

def plot_oer_energies(all_paths):
    """21가지 OER 경로의 반응 에너지를 비교하는 그래프"""
    if not all_paths:
        print("No valid paths found for plotting")
        return
    
    # 첫 번째 그래프: 각 경로별 막대 그래프
    fig1, ax1 = plt.subplots(1, 1, figsize=(20, 8))
    
    steps = ['Step 1', 'Step 2', 'Step 3', 'Step 4']
    colors = ['blue', 'red', 'green', 'magenta', 'orange', 'purple', 'brown', 'pink', 'gray', 
              'cyan', 'yellow', 'lime', 'navy', 'olive', 'teal', 'coral', 'gold', 'indigo', 
              'maroon', 'silver', 'tan']
    path_names = list(all_paths.keys())
    
    x = np.arange(len(steps))
    width = 0.03  # 21개 경로를 위해 막대 폭을 더 줄임
    
    for i, (path_name, energies) in enumerate(all_paths.items()):
        offset = (i - len(path_names)/2 + 0.3) * width
        bars = ax1.bar(x + offset, energies, width, label=path_name, color=colors[i % len(colors)], alpha=0.8, edgecolor='black', linewidth=0.3)
        
        # 각 막대 위에 값 표시 (모든 막대에 표시, 90도 회전)
        for j, (bar, energy) in enumerate(zip(bars, energies)):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                    f'{energy:.2f}', ha='center', va='bottom', fontsize=7, rotation=90)
    
    ax1.set_xlabel('Reaction Steps', fontsize=12)
    ax1.set_ylabel('Reaction Energy (eV)', fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(steps)
    ax1.axhline(y=1.23, color='black', linestyle='--', linewidth=1)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend(fontsize=10, ncol=6)
    
    plt.tight_layout()
    plt.savefig('OER_energies_comparison_4steps.png', dpi=300, bbox_inches='tight')
    # plt.show()
    
    # 두 번째 그래프: 각 경로의 최대 에너지 비교
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 8))
    
    max_energies = [max(energies) for energies in all_paths.values()]
    
    bars2 = ax2.bar(path_names, max_energies, color=colors[:len(path_names)], alpha=0.8, edgecolor='black', linewidth=1)
    ax2.axhline(y=1.23, color='black', linestyle='--', linewidth=1)
    
    # 각 막대 위에 값 표시
    for bar, energy in zip(bars2, max_energies):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                f'{energy:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax2.set_ylabel('Maximum Step Energy (eV)', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend(fontsize=10)
    
    # x축 라벨 회전
    ax2.tick_params(axis='x', rotation=90)
    
    plt.tight_layout()
    plt.savefig('OER_overpotential_comparison_4steps.png', dpi=300, bbox_inches='tight')
    # plt.show()

def plot_individual_energetics_diagrams(all_paths, oer_data):
    """각 path에 대해 개별 energetics diagram을 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting individual diagrams")
        return
    
    # 각 path별로 개별 그래프 생성
    for i, (path_name, energies) in enumerate(all_paths.items()):
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        
        # 데이터에서 해당 path 정보 찾기 (인덱스로 찾기)
        path_data = oer_data[i] if i < len(oer_data) else None
        if not path_data:
            continue
            
        # x축: 반응 단계 (1, 2, 3, 4)
        steps = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5]

        energy0 = 0.0
        energy1 = path_data['step1']
        energy2 = path_data['step2'] + energy1
        energy3 = path_data['step3'] + energy2
        energy4 = path_data['step4'] + energy3
        surface = path_data['surface']
        
        # y축: step energies
        step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4]
        
        # 선 그래프로 그리기
        ax.plot(steps, step_energies, '-', color='black', linewidth=2)

        # 그래프 설정
        ax.set_xticks([])
        ax.set_xlabel('Reaction Coordinate', fontsize=12)
        ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
        
        # y축 범위 설정 (최소값과 최대값에 여유 공간 추가)
        y_min = min(step_energies) - 0.2
        y_max = max(step_energies) + 0.3
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(0, 5)
                
        # 최대 에너지와 과전위 정보 추가
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'])
        # overpotential = int((max_energy - 1.23) * 1000)
        overpotential = max_energy - 1.23
        # ax.text(0.02, 0.97, f'Surface: {surface}\nΔG1: {path_data["step1"]:.2f} eV\nΔG2: {path_data["step2"]:.2f} eV\nΔG3: {path_data["step3"]:.2f} eV\nΔG4: {path_data["step4"]:.2f} eV\nOverpotential: {overpotential} mV', 
        ax.text(0.02, 0.97, f'Surface: {surface}\nΔG1: {path_data["step1"]:.2f} eV\nΔG2: {path_data["step2"]:.2f} eV\nΔG3: {path_data["step3"]:.2f} eV\nΔG4: {path_data["step4"]:.2f} eV\nOverpotential: {overpotential:.2f} V', 
                transform=ax.transAxes, fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        # 파일명에 path 정보 포함하여 저장
        filename = f'energetics_diagram_4steps_{path_name.lower()}_{path_data["surface"].lower()}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        # plt.show()
        plt.close()
        
        print(f"Saved: {filename}")

def plot_all_energetics_diagrams_combined(all_paths, oer_data):
    """모든 path의 energetics diagram을 하나의 그래프에 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting combined diagram")
        return
    
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    
    # 색상 팔레트
    colors = ['blue', 'red', 'green', 'magenta', 'orange', 'purple', 'brown', 'pink', 'gray', 
              'cyan', 'yellow', 'lime', 'navy', 'olive', 'teal', 'coral', 'gold', 'indigo', 
              'maroon', 'silver', 'tan']
    
    steps = [1, 2, 3, 4]
    
    for i, (path_name, energies) in enumerate(all_paths.items()):
        # 데이터에서 해당 path 정보 찾기 (인덱스로 찾기)
        path_data = oer_data[i] if i < len(oer_data) else None
        if not path_data:
            continue
            
        # step energies
        step_energies = [path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4']]
        
        # 선 그래프로 그리기
        ax.plot(steps, step_energies, 'o-', linewidth=1.5, markersize=4, 
                color=colors[i % len(colors)], alpha=0.8, 
                label=f'{path_name}: {path_data["surface"]}')
    
    # 열역학적 한계선 표시
    ax.axhline(y=1.23, color='red', linestyle='--', linewidth=2, alpha=0.7, 
              label='Thermodynamic limit (1.23 V)')
    
    # 그래프 설정
    ax.set_xlabel('Reaction Steps', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reaction Energy (eV)', fontsize=12, fontweight='bold')
    ax.set_title('OER Energetics Diagrams - All Paths Combined', fontsize=16, fontweight='bold')
    
    # x축 설정
    ax.set_xticks(steps)
    ax.set_xticklabels([f'Step {i}' for i in steps])
    
    # 범례 (2열로 배치)
    ax.legend(fontsize=8, ncol=2, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    
    plt.tight_layout()
    plt.savefig('all_energetics_diagrams_combined_4steps.png', dpi=300, bbox_inches='tight')
    # plt.show()
    
    print("Saved: all_energetics_diagrams_combined_4steps.png")
    
    # 결과를 텍스트 파일로 저장
    with open('OER_energies_comparison_4steps.txt', 'w') as f:
        f.write("OER Energetics Comparison - All Paths\n")
        f.write("=" * 60 + "\n\n")
        
        for path_name, energies in all_paths.items():
            f.write(f"{path_name}:\n")
            f.write(f"  Step 1: {energies[0]:.6f} eV\n")
            f.write(f"  Step 2: {energies[1]:.6f} eV\n")
            f.write(f"  Step 3: {energies[2]:.6f} eV\n")
            f.write(f"  Step 4: {energies[3]:.6f} eV\n")
            f.write(f"  Maximum: {max(energies):.6f} eV\n")
            f.write(f"  Overpotential: {max(energies) - 1.23:.6f} eV\n\n")
        
        f.write("Summary:\n")
        f.write("-" * 30 + "\n")
        for path_name, energies in all_paths.items():
            f.write(f"{path_name}: Max = {max(energies):.3f} eV, Overpotential = {max(energies) - 1.23:.3f} eV\n")
        
        # 최적 경로 찾기
        best_path = min(all_paths.keys(), key=lambda x: max(all_paths[x]))
        f.write(f"\nBest path: {best_path} (lowest overpotential)\n")

if __name__ == "__main__":
    all_paths, oer_data = calculate_oer_energies()
    
    if not all_paths:
        print("OER 에너지 계산에 실패했습니다. 경로와 파일을 확인하세요.")
    else:
        # 표 형태로 데이터 출력
        print_oer_table(oer_data)
        
        # 막대 그래프 출력
        plot_oer_energies(all_paths)
        
        # 개별 energetics diagram 출력
        print("\nGenerating individual energetics diagrams...")
        plot_individual_energetics_diagrams(all_paths, oer_data)
        
        # 모든 energetics diagram을 하나의 그래프에 출력
        print("\nGenerating combined energetics diagram...")
        plot_all_energetics_diagrams_combined(all_paths, oer_data)