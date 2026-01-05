from ase.io import read
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
import numpy as np
import argparse

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

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
gh2o, gh2 = -14.586115, -6.905510 

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

t = 298.15
dgh = (1.577213-t*0.000460)/8
dgoh = (2.928444-t*0.001666)/8
dgo = (0.690383-t*0.000800)/8

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

def load_label_csv(label_csv_path):
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
            print(f"Loaded OH counts from label.csv")
        except Exception as e:
            print(f"Warning: Could not read label.csv: {e}")
    else:
        print(f"Warning: label.csv not found. OH counts will be set to 0.")
    return file_oh_counts

def calculate_gibbs_correction(n_H, n_O, oh_count):
    """Gibbs free energy 보정을 계산하는 함수"""
    correction = (
        - (n_H - oh_count) * (gh - dgh)  # H atoms (excluding those in OH)
        - (n_O - oh_count) * (go - dgo)   # O atoms (excluding those in OH)
        - oh_count * (goh - dgoh)         # OH groups
    )
    return correction

def calculate_charge(n_H, n_O, oh_count):
    """구조의 전하를 계산하는 함수"""
    charge = (n_H - oh_count) * 1 + (n_O - oh_count) * (-2) + oh_count * (-1)
    return charge

def calculate_electrons_needed(prev_n_H, prev_n_O, prev_oh_count, 
                                curr_n_H, curr_n_O, curr_oh_count):
    """반응에 필요한 전자 수를 계산하는 함수"""
    prev_charge = calculate_charge(prev_n_H, prev_n_O, prev_oh_count)
    curr_charge = calculate_charge(curr_n_H, curr_n_O, curr_oh_count)
    charge_change = curr_charge - prev_charge
    n_electrons = -charge_change
    return n_electrons, charge_change

def process_path(path, file_order, voltage=0.0):
    """주어진 경로에서 데이터를 처리하고 그래프 데이터를 반환하는 함수"""
    # label.csv 파일 경로
    label_csv_path = os.path.join(path, 'label.csv')
    file_oh_counts = load_label_csv(label_csv_path)
    
    # 에너지 읽기 및 보정 적용
    energies = []
    corrected_energies = []
    reaction_energies = []
    cumulative_energies = []
    labels = []
    atom_counts_list = []
    
    for filename, label in file_order:
        filepath = os.path.join(path, filename)
        if os.path.exists(filepath):
            energy = get_energy_from_json(filepath)
            if energy is not None:
                n_H, n_O = count_atoms_from_json(filepath)
                oh_count = file_oh_counts.get(filename, 0)
                gibbs_correction = calculate_gibbs_correction(n_H, n_O, oh_count)
                corrected_energy = energy + gibbs_correction
                
                energies.append(energy)
                corrected_energies.append(corrected_energy)
                labels.append(label)
                atom_counts_list.append((n_H, n_O, oh_count))
            else:
                print(f"Warning: Could not read energy from {filepath}")
    
    if not corrected_energies:
        print(f"No energy data found in {path}!")
        return None
    
    # 각 단계의 반응 에너지 계산
    n_sites = 8  # 활성사이트 개수
    reaction_potentials = []
    
    for i in range(len(corrected_energies)):
        if i == 0:
            reaction_energy = 0.0
            cumulative_energy = 0.0
            reaction_potential = None
            print(f"Step 0: {labels[0]} (기준점, 에너지 = 0.00 eV/site)")
        else:
            prev_n_H, prev_n_O, prev_oh_count = atom_counts_list[i-1]
            curr_n_H, curr_n_O, curr_oh_count = atom_counts_list[i]
            
            n_electrons_total, charge_change = calculate_electrons_needed(
                prev_n_H, prev_n_O, prev_oh_count,
                curr_n_H, curr_n_O, curr_oh_count
            )
            
            n_electrons = n_electrons_total / n_sites
            reaction_energy_0V = (corrected_energies[i] - corrected_energies[i-1]) / n_sites
            
            if abs(charge_change) > 0.01:
                reaction_energy = reaction_energy_0V - n_electrons * voltage
                reaction_potential = reaction_energy_0V / n_electrons if n_electrons > 0 else None
            else:
                reaction_energy = reaction_energy_0V
                reaction_potential = None
            
            cumulative_energy = cumulative_energies[i-1] + reaction_energy
            
            print(f"Step {i}: {labels[i-1]} → {labels[i]}")
            print(f" 반응 에너지 (0V 기준) = {reaction_energy_0V:.2f} eV/site")
        
        reaction_energies.append(reaction_energy)
        cumulative_energies.append(cumulative_energy)
        reaction_potentials.append(reaction_potential)
    
    # labels 확장
    expanded_labels = []
    for label in labels:
        expanded_labels.append(label)
    labels = expanded_labels
    
    relative_energies = cumulative_energies
    relative_energies[0] = 0.0
    
    return {
        'relative_energies': relative_energies,
        'labels': labels,
        'reaction_energies': reaction_energies,
        'reaction_potentials': reaction_potentials,
        'path': path
    }

def plot_energetics_dual(path1, path2, file_order, voltage=0.0, ymin=None, ymax=None, 
                         xmin=None, xmax=None, show_legend=False, show_text=True,
                         label1=None, label2=None, color1='blue', color2='red',
                         zorder1=1, zorder2=2, suffix=None, figsize=(6, 4), gap=0.5):
    """두 경로의 energetics 다이어그램을 한 이미지에 그리는 함수
    
    Parameters:
    -----------
    path1 : str
        첫 번째 경로
    path2 : str
        두 번째 경로
    file_order : list
        파일 이름과 표시 이름 매핑 리스트
    voltage : float, optional
        작동전압 (V). 기본값은 0.0 V.
    ymin : float, optional
        y축 최소값
    ymax : float, optional
        y축 최대값
    xmin : float, optional
        x축 최소값
    xmax : float, optional
        x축 최대값
    show_legend : bool, optional
        범례 표시 여부
    show_text : bool, optional
        각 점의 레이블 표시 여부
    label1 : str, optional
        첫 번째 경로의 레이블 (범례용)
    label2 : str, optional
        두 번째 경로의 레이블 (범례용)
    color1 : str, optional
        첫 번째 경로의 색상
    color2 : str, optional
        두 번째 경로의 색상
    zorder1 : int, optional
        첫 번째 경로의 zorder (그리는 순서, 높을수록 위에 표시)
    zorder2 : int, optional
        두 번째 경로의 zorder (그리는 순서, 높을수록 위에 표시)
    suffix : str, optional
        파일 이름 뒤에 붙일 접미사
    figsize : tuple, optional
        그래프 크기 (width, height) 인치 단위. 기본값은 (6, 4)
    """
    # 두 경로에서 데이터 처리
    print(f"\n=== Processing Path 1 ===")
    data1 = process_path(path1, file_order, voltage)
    print(f"\n=== Processing Path 2 ===")
    data2 = process_path(path2, file_order, voltage)
    
    if data1 is None or data2 is None:
        print("Error: Could not process one or both paths!")
        return
    
    # 그래프 그리기
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # 두 경로의 데이터
    datasets = [
        (data1, color1, label1, zorder1),
        (data2, color2, label2, zorder2)
    ]
    
    # y축 범위 계산 (두 경로 모두 포함)
    all_energies = data1['relative_energies'] + data2['relative_energies']
    if ymin is None:
        y_min = min(all_energies) - 0.2
    else:
        y_min = ymin
    
    if ymax is None:
        y_max = max(all_energies) + 0.3
    else:
        y_max = ymax
    
    # 각 경로에 대해 그래프 그리기
    for data, color, path_label, zorder in datasets:
        relative_energies = data['relative_energies']
        labels = data['labels']
        reaction_energies = data['reaction_energies']
        reaction_potentials = data['reaction_potentials']
        
        x_positions = list(range(len(relative_energies)))
        
        # 선 그리기
        electrochemical_plotted = False
        chemical_plotted = False
        
        for i in range(len(x_positions)):
            x_segment1 = [x_positions[i], x_positions[i]+gap]
            y_segment1 = [relative_energies[i], relative_energies[i]]
            ax.plot(x_segment1, y_segment1, linestyle='-', color=color, 
                   linewidth=2, zorder=zorder)

            if i != len(x_positions) - 1:
                x_segment2 = [x_positions[i]+gap, x_positions[i+1]]
                y_segment2 = [relative_energies[i], relative_energies[i+1]]
                if reaction_potentials[i+1] is not None and reaction_potentials[i+1] != 0:
                    linestyle = '-'
                    label = f'{path_label} (electrochemical step)' if not electrochemical_plotted else ''
                    electrochemical_plotted = True
                else:
                    linestyle = '--'
                    label = f'{path_label} (chemical step)' if not chemical_plotted else ''
                    chemical_plotted = True
                ax.plot(x_segment2, y_segment2, linestyle=linestyle, color=color, 
                    label=label, linewidth=1, zorder=zorder)
        
        # 텍스트 레이블 추가
        if show_text:
            for i, (x, y, label) in enumerate(zip(x_positions, relative_energies, labels)):
                if i < len(reaction_energies) - 1:
                    next_reaction_energy = reaction_energies[i+1]
                    if next_reaction_energy > 0.1:
                        y_offset = -0.05
                        ha = 'right'
                        va = 'top'
                    else:
                        y_offset = 0.05
                        ha = 'left'
                        va = 'bottom'
                else:
                    y_offset = 0.05
                    ha = 'left'
                    va = 'bottom'
                
                ax.text(x, y + y_offset, label, ha=ha, va=va, fontsize=8, 
                       rotation=45, color=color, alpha=0.7)
    
    # 범례 추가
    if show_legend:
        ax.legend(loc='best')
    
    # 그래프 설정
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.xaxis.set_label_position('top')  # xlabel을 위에 표시
    ax.set_ylabel('Relative Energy (ΔG, eV/site)', fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    
    # x축 범위 설정
    max_len = max(len(data1['relative_energies']), len(data2['relative_energies']))
    if xmin is None:
        x_min = 0
    else:
        x_min = xmin
    
    if xmax is None:
        x_max = max_len - 1
    else:
        x_max = xmax
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    if suffix:
        output_filename = f'energetics_diagram_dual_{suffix}.png'
    else:
        output_filename = 'energetics_diagram_dual.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\nDual energetics diagram saved as '{output_filename}'")
    # print(f"Note: Energies are corrected with Gibbs free energy corrections")
    # print(f"Note: Energies are normalized per active site (n_sites = 8)")
    
    plt.show()

if __name__ == '__main__':
    # 파일 이름과 표시 이름 매핑
    file_order = [
        ('Ir_Ir.json', 'Ir-Ir'),
        ('OH_Ir_Ir.json', 'OH-Ir-Ir'),
        ('O_Ir_Ir.json', 'O-Ir-Ir'),
        ('Ir_O_Ir.json', 'Ir-O-Ir'),
        ('OH_Ir_O_Ir.json', 'OH-Ir-O-Ir'),
        ('O_Ir_O_Ir.json', 'O-Ir-O-Ir'),
        ('ir_O_ir_O_Ir.json', '½Ir-O-½Ir-O-Ir'),
        ('Ir_O_O_Ir.json', 'Ir-O-O-Ir'),
        ('OH_Ir_O_O_Ir.json', 'OH-Ir-O-O-Ir'),
        ('O_Ir_O_O_Ir.json', 'O-Ir-O-O-Ir'),
        ('O_ir_O_ir_O_Ir.json', 'O-½Ir-O-½Ir-O-Ir'),
    ]
    
    # 명령줄 인자 파싱
    parser = argparse.ArgumentParser(description='Plot dual energetics diagram from two paths')
    parser.add_argument('path1', type=str, nargs='?', 
                        default='/Users/jiuy97/Desktop/irfe/revision/1_Ir_Ir_Ir/pourbaix',
                        help='첫 번째 경로 (기본값: /Users/jiuy97/Desktop/irfe/revision/1_Ir_Ir_Ir/pourbaix)')
    parser.add_argument('path2', type=str, nargs='?',
                        default='/Users/jiuy97/Desktop/irfe/revision/2_Fe_Ir_Ir/pourbaix',
                        help='두 번째 경로 (기본값: /Users/jiuy97/Desktop/irfe/revision/2_Fe_Ir_Ir/pourbaix)')
    parser.add_argument('--voltage', '-V', type=float, default=0.0,
                        help='작동전압 (V). 기본값: 0.0 V')
    parser.add_argument('--ymin', type=float, default=None,
                        help='y축 최소값')
    parser.add_argument('--ymax', type=float, default=None,
                        help='y축 최대값')
    parser.add_argument('--xmin', type=float, default=None,
                        help='x축 최소값')
    parser.add_argument('--xmax', type=float, default=None,
                        help='x축 최대값')
    parser.add_argument('--show-legend', action='store_true',
                        help='범례를 표시합니다.')
    parser.add_argument('--show-text', dest='show_text', action='store_const', const=True, default=True,
                        help='각 점의 레이블을 표시합니다. (기본값: True)')
    parser.add_argument('--no-show-text', dest='show_text', action='store_const', const=False,
                        help='각 점의 레이블을 표시하지 않습니다.')
    parser.add_argument('--label1', type=str, default='Ir',
                        help='첫 번째 경로의 레이블 (범례용)')
    parser.add_argument('--label2', type=str, default='IrFe',
                        help='두 번째 경로의 레이블 (범례용)')
    parser.add_argument('--color1', type=str, default='black',
                        help='첫 번째 경로의 색상 (기본값: black)')
    parser.add_argument('--color2', type=str, default='red',
                        help='두 번째 경로의 색상 (기본값: red)')
    parser.add_argument('--zorder1', type=int, default=2,
                        help='첫 번째 경로의 zorder (그리는 순서, 높을수록 위에 표시, 기본값: 2)')
    parser.add_argument('--zorder2', type=int, default=1,
                        help='두 번째 경로의 zorder (그리는 순서, 높을수록 위에 표시, 기본값: 1)')
    parser.add_argument('--suffix', type=str, default=None,
                        help='파일 이름 뒤에 붙일 접미사 (예: --suffix test → energetics_diagram_dual_test.png)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[6, 4],
                        metavar=('WIDTH', 'HEIGHT'),
                        help='그래프 크기 (width height) 인치 단위. 기본값: 6 4')
    parser.add_argument('--gap', type=float, default=0.5)
    
    args = parser.parse_args()
    
    plot_energetics_dual(
        path1=args.path1,
        path2=args.path2,
        file_order=file_order,
        voltage=args.voltage,
        ymin=args.ymin,
        ymax=args.ymax,
        xmin=args.xmin,
        xmax=args.xmax,
        show_legend=args.show_legend,
        show_text=args.show_text,
        label1=args.label1,
        label2=args.label2,
        color1=args.color1,
        color2=args.color2,
        zorder1=args.zorder1,
        zorder2=args.zorder2,
        suffix=args.suffix,
        figsize=tuple(args.figsize),
        gap=args.gap
    )
