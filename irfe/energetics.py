from ase.io import read
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
import numpy as np
import argparse

# ========================================
# THERMODYNAMIC CONSTANTS (from HybridPB/pourbaix.py)
# ========================================

# Gas phase reference energies (DFT calculated)
h2 = -6.77149190   # H₂ molecule total energy (eV)
h2o = -14.23091949 # H₂O molecule total energy (eV)

# Thermodynamic corrections for gas phases
zpeh2o = 0.558  # Zero-point energy correction (eV)
cvh2o = 0.103   # Heat capacity correction (eV)
tsh2o = 0.675   # Entropy correction T×S (eV)

zpeh2 = 0.268   # Zero-point energy correction (eV)
cvh2 = 0.0905   # Heat capacity correction (eV)
tsh2 = 0.408    # Entropy correction T×S (eV)

# Gibbs free energies of gas phases
gh2o = h2o + zpeh2o - tsh2o + cvh2o  # G(H₂O) gas phase
gh2 = h2 + zpeh2 - tsh2 + cvh2       # G(H₂) gas phase

# Derived chemical potentials
gh = gh2 / 2                    # μ(H) = ½G(H₂)
go = gh2o - gh2                 # μ(O) = G(H₂O) - G(H₂)
goh = gh2o - gh2 / 2            # μ(OH) = G(H₂O) - ½G(H₂)

# Adsorbate thermodynamic corrections
zpeoh = 0.376   # Zero-point energy (eV)
cvoh = 0.042    # Heat capacity (eV)
tsoh = 0.066    # Entropy T×S (eV)

zpeo = 0.064    # Zero-point energy (eV)
cvo = 0.034     # Heat capacity (eV)
tso = 0.060     # Entropy T×S (eV)

# Adsorbate Gibbs free energy corrections
dgo = zpeo + cvo - tso        # Gibbs correction for O* adsorbate
dgoh = zpeoh + cvoh - tsoh    # Gibbs correction for OH* adsorbate
dgh = dgoh - dgo              # Gibbs correction for H* adsorbate (derived)

def get_energy_from_json(json_path):
    """JSON 파일에서 에너지를 읽어오는 함수"""
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def count_atoms_from_json(json_path):
    """JSON 파일에서 H, O 원자 개수를 세는 함수"""
    try:
        atoms = read(json_path)
        symbols = atoms.get_chemical_symbols()
        n_H = symbols.count('H')
        n_O = symbols.count('O')
        return n_H, n_O
    except Exception as e:
        print(f"Error counting atoms from {json_path}: {e}")
        return 0, 0

def load_label_csv(label_csv_path='label.csv'):
    """label.csv 파일을 읽어서 OH 개수 정보를 반환하는 함수"""
    file_oh_counts = {}
    if os.path.exists(label_csv_path):
        try:
            # CSV 파일 읽기: json_name, label, #OH 컬럼
            label_df = pd.read_csv(label_csv_path, header=None, 
                                  names=['json_name', 'label', '#OH'])
            for idx, row in label_df.iterrows():
                json_name = row['json_name']
                if '#OH' in row and not pd.isna(row['#OH']):
                    file_oh_counts[json_name] = float(row['#OH'])
            print(f"Loaded OH counts from {label_csv_path}")
        except Exception as e:
            print(f"Warning: Could not read {label_csv_path}: {e}")
    else:
        print(f"Warning: {label_csv_path} not found. OH counts will be set to 0.")
    return file_oh_counts

def calculate_gibbs_correction(n_H, n_O, oh_count):
    """Gibbs free energy 보정을 계산하는 함수"""
    # formation_energy_correction = -(H - oh_count) * (gh - dgh) - (O - oh_count) * (go - dgo) - oh_count * (goh - dgoh)
    correction = (
        - (n_H - oh_count) * (gh - dgh)  # H atoms (excluding those in OH)
        - (n_O - oh_count) * (go - dgo)   # O atoms (excluding those in OH)
        - oh_count * (goh - dgoh)         # OH groups
    )
    return correction

def calculate_charge(n_H, n_O, oh_count):
    """구조의 전하를 계산하는 함수"""
    # H: +1, O: -2, OH: -1 (H+O = +1-2 = -1)
    # OH에 포함된 H와 O는 이미 OH로 계산되므로 제외
    charge = (n_H - oh_count) * 1 + (n_O - oh_count) * (-2) + oh_count * (-1)
    return charge

def calculate_electrons_needed(prev_n_H, prev_n_O, prev_oh_count, 
                                curr_n_H, curr_n_O, curr_oh_count):
    """반응에 필요한 전자 수를 계산하는 함수"""
    prev_charge = calculate_charge(prev_n_H, prev_n_O, prev_oh_count)
    curr_charge = calculate_charge(curr_n_H, curr_n_O, curr_oh_count)
    
    # 전하 변화 (양수면 전자가 필요, 음수면 전자가 방출)
    charge_change = curr_charge - prev_charge
    
    # 전자 수 = -전하 변화 (전자가 필요하면 양수)
    n_electrons = -charge_change
    
    return n_electrons, charge_change


def plot_energetics_diagram(voltage=0.0, ymin=None, ymax=None, xmin=None, xmax=None):
    """Energetics 다이어그램을 그리는 함수
    
    Parameters:
    -----------
    voltage : float, optional
        작동전압 (V). 기본값은 0.0 V.
        전기화학 반응의 반응 에너지에 반영됩니다.
    ymin : float, optional
        y축 최소값. 기본값은 None (자동 계산: min(relative_energies) - 0.2)
    ymax : float, optional
        y축 최대값. 기본값은 None (자동 계산: max(relative_energies) + 0.3)
    xmin : float, optional
        x축 최소값. 기본값은 None (자동 계산: 0)
    xmax : float, optional
        x축 최대값. 기본값은 None (자동 계산: len(x_positions) - 1)
    """
    # label.csv 파일에서 OH 개수 읽기
    file_oh_counts = load_label_csv('label.csv')
    
    # 파일 이름과 표시 이름 매핑
    file_order = [
        ('Ir_Ir.json', 'Ir-Ir'),
        ('OH_Ir_Ir.json', 'OH-Ir-Ir'),
        ('O_Ir_Ir.json', 'O-Ir-Ir'),
        ('Ir_O_Ir.json', 'Ir-O-Ir'),
        ('OH_Ir_O_Ir.json', 'OH-Ir-O-Ir'),
        ('O_Ir_O_Ir.json', 'O-Ir-O-Ir'),
        ('Ir_O_O_Ir.json', 'Ir-O-O-Ir'),
        ('OH_Ir_O_O_Ir.json', 'OH-Ir-O-O-Ir'),
        ('O_Ir_O_O_Ir.json', 'O-Ir-O-O-Ir'),
    ]
    
    # 에너지 읽기 및 보정 적용
    energies = []
    corrected_energies = []
    reaction_energies = []  # 각 단계의 반응 에너지
    cumulative_energies = []  # 누적 반응 에너지
    labels = []
    atom_counts_list = []  # 각 구조의 원자 개수 저장 (H, O, OH)
    
    for filename, label in file_order:
        if os.path.exists(filename):
            energy = get_energy_from_json(filename)
            if energy is not None:
                # 원자 개수 세기
                n_H, n_O = count_atoms_from_json(filename)
                
                # OH 개수 가져오기 (label.csv에서, 없으면 0)
                oh_count = file_oh_counts.get(filename, 0)
                
                # Gibbs free energy 보정 계산
                gibbs_correction = calculate_gibbs_correction(n_H, n_O, oh_count)
                
                # 보정된 에너지
                corrected_energy = energy + gibbs_correction
                
                energies.append(energy)
                corrected_energies.append(corrected_energy)
                labels.append(label)
                atom_counts_list.append((n_H, n_O, oh_count))
                
                print(f"{label}: E_DFT = {energy:.6f} eV, "
                      f"H={n_H}, O={n_O}, OH={oh_count}, "
                      f"G_corr = {gibbs_correction:.6f} eV, "
                      f"E_corrected = {corrected_energy:.6f} eV")
            else:
                print(f"Warning: Could not read energy from {filename}")
        else:
            print(f"Warning: File {filename} not found")
    
    if not corrected_energies:
        print("No energy data found!")
        return
    
    # 각 단계의 반응 에너지 계산 (보정된 에너지의 차이만 사용)
    # calculate_gibbs_correction에서 이미 gh, go, goh를 반영했으므로
    # 단순히 보정된 에너지의 차이만 계산하면 됨
    # 활성사이트가 셀당 8개이므로 반응 에너지를 8로 정규화
    n_sites = 8  # 활성사이트 개수
    
    print("\n=== 반응 에너지 및 반응 포텐셜 계산 (활성사이트당 정규화) ===")
    reaction_potentials = []  # 각 단계의 반응 포텐셜
    
    for i in range(len(corrected_energies)):
        if i == 0:
            # 첫 번째 단계는 기준점 (에너지 = 0)
            reaction_energy = 0.0
            cumulative_energy = 0.0
            n_electrons = 0
            reaction_potential = None
            print(f"Step 0: {labels[0]} (기준점, 에너지 = 0.0 eV/site)")
        else:
            # 이전 구조와 현재 구조의 원자 개수
            prev_n_H, prev_n_O, prev_oh_count = atom_counts_list[i-1]
            curr_n_H, curr_n_O, curr_oh_count = atom_counts_list[i]
            
            # 필요한 전자 수 계산
            n_electrons_total, charge_change = calculate_electrons_needed(
                prev_n_H, prev_n_O, prev_oh_count,
                curr_n_H, curr_n_O, curr_oh_count
            )
            
            # 활성사이트당 전자 수 정규화 (8로 나눔)
            n_electrons = n_electrons_total / n_sites
            
            # 보정된 에너지의 차이 = 반응 에너지 (0 V 기준)
            # 활성사이트당 정규화 (8로 나눔)
            reaction_energy_0V = (corrected_energies[i] - corrected_energies[i-1]) / n_sites
            
            # 전기화학 반응인지 확인 (전하 변화가 있으면 전기화학 반응)
            if abs(charge_change) > 0.01:  # 전하 변화가 있으면 전기화학 반응
                # 전기화학 반응: 작동전압 반영
                # ΔG = ΔG_0 - n_e * e * U (e = 1 eV/V)
                reaction_energy = reaction_energy_0V - n_electrons * voltage
                # 반응 포텐셜 = 반응 에너지 / 전자 수 (활성사이트당)
                reaction_potential = reaction_energy_0V / n_electrons if n_electrons > 0 else None
            else:
                # 화학 반응 (전하 변화 없음): 전압 영향 없음
                reaction_energy = reaction_energy_0V
                reaction_potential = None
            
            cumulative_energy = cumulative_energies[i-1] + reaction_energy
            
            print(f"Step {i}: {labels[i-1]} → {labels[i]}")
            # print(f"  필요한 전자 수 = {n_electrons:.2f} e⁻/site")
            print(f"  반응 에너지 (0V 기준) = {reaction_energy_0V:.6f} eV/site")
            # if abs(charge_change) > 0.01:
            #     print(f"  작동전압 반영: {reaction_energy:.6f} eV/site (U = {voltage:.3f} V)")
            # print(f"  누적 에너지 = {cumulative_energy:.6f} eV/site")
            # if reaction_potential is not None:
            #     print(f"  반응 포텐셜 = {reaction_potential:.6f} V")
            # else:
            #     print(f"  반응 포텐셜 = N/A (화학 반응)")
        
        reaction_energies.append(reaction_energy)
        cumulative_energies.append(cumulative_energy)
        reaction_potentials.append(reaction_potential)
    
    # 첫 번째 에너지를 기준으로 상대 에너지 계산 (첫 번째 단계 = 0)
    relative_energies = cumulative_energies
    
    # 명시적으로 첫 번째 단계를 0으로 설정 (확인용)
    relative_energies[0] = 0.0
    
    # 그래프 그리기
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    
    # x축: 반응 단계 (0부터 시작)
    x_positions = list(range(len(relative_energies)))
    
    # 전기화학 반응과 화학 반응을 선 스타일로 구분하여 그리기
    # 전기화학 반응: solid line (-)
    # 화학 반응: dashed line (--)
    # 모든 색상은 검정색
    line_color = 'black'
    
    # 범례를 위한 플래그
    electrochemical_plotted = False
    chemical_plotted = False
    
    # 각 선분을 개별적으로 그리기
    for i in range(len(x_positions) - 1):
        x_segment = [x_positions[i], x_positions[i+1]]
        y_segment = [relative_energies[i], relative_energies[i+1]]
        
        # 전기화학 반응인지 화학 반응인지 확인
        if reaction_potentials[i+1] is not None:
            # 전기화학 반응: solid line
            linestyle = '-'
            label = 'Electrochemical' if not electrochemical_plotted else ''
            electrochemical_plotted = True
        else:
            # 화학 반응: dashed line
            linestyle = '--'
            label = 'Chemical' if not chemical_plotted else ''
            chemical_plotted = True
        
        ax.plot(x_segment, y_segment, linestyle=linestyle, color=line_color, label=label)
    
    # 모든 점에 동일한 검정색 마커 추가
    for i, (x, y) in enumerate(zip(x_positions, relative_energies)):
        ax.plot(x, y, 'o', color='black', markerfacecolor='white', markeredgewidth=1.5, zorder=5)
    
    # 각 점에 레이블 추가 (다음 단계의 반응 에너지에 따라 위치 조정)
    for i, (x, y, label) in enumerate(zip(x_positions, relative_energies, labels)):
        # 다음 단계의 반응 에너지 확인
        if i < len(reaction_energies) - 1:
            # 다음 단계(i+1)의 반응 에너지 확인
            next_reaction_energy = reaction_energies[i+1]
            
            if next_reaction_energy > 0:
                # 다음 단계의 반응 에너지가 양수: 아래쪽에 표시
                y_offset = -0.05
                ha = 'right'
                va = 'top'
            else:
                # 다음 단계의 반응 에너지가 0 이하: 위쪽에 표시
                y_offset = 0.05
                ha = 'left'
                va = 'bottom'
        else:
            # 마지막 점: 위쪽에 표시
            y_offset = 0.05
            ha = 'left'
            va = 'bottom'
        
        ax.text(x, y + y_offset, label, ha=ha, va=va, fontsize=9, rotation=45)
    
    # 범례 추가
    # ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0))
    ax.legend(loc='best')
    
    # 그래프 설정
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.set_ylabel('Relative Gibbs Free Energy (eV/site)', fontsize=12)
    # ax.set_title('Energetics Diagram (with Gibbs Free Energy Correction, normalized per site)', fontsize=14, fontweight='bold')
    # ax.grid(True, alpha=0.3, linestyle='--')
    
    # x축 레이블 설정
    ax.set_xticks([])
    # ax.set_xticks(x_positions)
    # ax.set_xticklabels([f'Step {i+1}' for i in range(len(x_positions))], rotation=0)
    
    # x축 범위 설정
    if xmin is None:
        x_min = 0
    else:
        x_min = xmin
    
    if xmax is None:
        x_max = len(x_positions) - 1
    else:
        x_max = xmax
    
    ax.set_xlim(x_min, x_max)
    
    # y축 범위 설정 (최소값과 최대값에 여유 공간 추가)
    if ymin is None:
        y_min = min(relative_energies) - 0.2
    else:
        y_min = ymin
    
    if ymax is None:
        y_max = max(relative_energies) + 0.3
    else:
        y_max = ymax
    
    ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    plt.savefig('energetics_diagram.png', dpi=300, bbox_inches='tight')
    print(f"\nEnergetics diagram saved as 'energetics_diagram.png'")
    print(f"Note: Energies are corrected with Gibbs free energy corrections from label.csv")
    print(f"Note: Energies are normalized per active site (n_sites = {n_sites})")
    
    # 반응 포텐셜 요약 출력
    print("\n=== 반응 포텐셜 요약 ===")
    for i in range(1, len(reaction_potentials)):
        if reaction_potentials[i] is not None:
            print(f"Step {i} ({labels[i-1]} → {labels[i]}): {reaction_potentials[i]:.4f} V")
    
    plt.show()

if __name__ == '__main__':
    # 명령줄 인자 파싱
    parser = argparse.ArgumentParser(description='Plot energetics diagram with Gibbs free energy correction')
    parser.add_argument('--voltage', '-V', type=float, default=0.0,
                        help='작동전압 (V). 전기화학 반응의 반응 에너지에 반영됩니다. (기본값: 0.0 V)')
    parser.add_argument('--ymin', type=float, default=None,
                        help='y축 최소값. 기본값은 자동 계산 (min(relative_energies) - 0.2)')
    parser.add_argument('--ymax', type=float, default=None,
                        help='y축 최대값. 기본값은 자동 계산 (max(relative_energies) + 0.3)')
    parser.add_argument('--xmin', type=float, default=None,
                        help='x축 최소값. 기본값은 자동 계산 (0)')
    parser.add_argument('--xmax', type=float, default=None,
                        help='x축 최대값. 기본값은 자동 계산 (len(x_positions) - 1)')
    
    args = parser.parse_args()
    
    plot_energetics_diagram(voltage=args.voltage, ymin=args.ymin, ymax=args.ymax, 
                           xmin=args.xmin, xmax=args.xmax)
