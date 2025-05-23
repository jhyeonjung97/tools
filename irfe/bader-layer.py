import os
import json
import glob
from ase.io import read
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.geometry import get_distances
from collections import defaultdict

# Base directory
base_dir = '/Users/hailey/Desktop/4_IrFe3'

# Initialize results list
results = []

# 사용할 표면 리스트
selected_systems = ['5_IrMn', '6_IrFe', '7_IrCo', '8_IrNi', '0_Ir']
system_names = ['IrMn', 'IrFe', 'IrCo', 'IrNi', 'Ir']

# 사용할 흡착물 리스트
selected_adsorbates = ['1_H', '2_OH', '3_O']
coverages = ['layer']

# Find all atoms_bader_charge.json files and sort them (선택된 표면, 흡착물만)
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = [f for f in sorted(glob.glob(pattern))
              if any(f'/{sys}/' in f for sys in selected_systems)
              and any(f'/{ads}/' in f for ads in selected_adsorbates)]

# Find slab atoms_bader_charge.json files (선택된 표면만)
slab_pattern = os.path.join(base_dir, 'slab/*_*/atoms_bader_charge.json')
slab_json_files = [f for f in sorted(glob.glob(slab_pattern)) if any(f'/{sys}/' in f for sys in selected_systems)]

def get_layer_charges_and_distances(atoms):
    # Get z-coordinates and sort atoms
    z_coords = atoms.get_positions()[:, 2]
    sorted_indices = np.argsort(z_coords)
    sorted_atoms = atoms[sorted_indices]
    
    # 고정된 크기(8개씩)로 레이어 분할
    layer_size = 8
    num_atoms = len(atoms)
    # 레이어 수 계산 (원자 수가 8의 배수가 아닐 경우 마지막 레이어는 더 적은 원자를 가질 수 있음)
    num_layers = int(np.ceil(num_atoms / layer_size))
    
    # 각 레이어에 속하는 원자 인덱스 계산
    layer_indices = []
    for i in range(num_layers):
        start_idx = i * layer_size
        end_idx = min((i + 1) * layer_size, num_atoms)
        layer_indices.append(np.arange(start_idx, end_idx))
    
    # Initialize layer charges and distances
    layer_charges = {}
    layer_distances = {}
    layer_z_avg = {}  # 각 레이어의 평균 z좌표 저장
    
    # Process each layer
    for layer_idx, indices in enumerate(layer_indices):
        # 인덱스가 1부터 시작하도록 layer 번호 조정 (L1, L2, L3, L4...)
        layer_num = layer_idx + 1
        
        layer_atoms = sorted_atoms[indices]
        layer_z_positions = layer_atoms.get_positions()[:, 2]
        layer_z_avg[f'L{layer_num}_z'] = round(np.mean(layer_z_positions), 2)  # 레이어의 평균 z좌표 저장
        
        # Calculate average Bader charges for Ir and TM
        ir_charges = []
        tm_charges = []
        ir_indices = []
        tm_indices = []
        
        for i, atom in enumerate(layer_atoms):
            symbol = atom.symbol
            charge = atom.charge
            
            if symbol == 'Ir':
                ir_charges.append(charge)
                ir_indices.append(i)
            elif symbol in ['Mn', 'Fe', 'Co', 'Ni']:
                tm_charges.append(charge)
                tm_indices.append(i)
        
        # Store average charges for this layer
        layer_charges[f'L{layer_num}_Ir'] = round(sum(ir_charges) / len(ir_charges), 2) if ir_charges else None
        layer_charges[f'L{layer_num}_TM'] = round(sum(tm_charges) / len(tm_charges), 2) if tm_charges else None
        
        # Calculate distances within the layer
        if ir_indices and tm_indices:
            # Calculate Ir-TM distances
            ir_positions = layer_atoms[ir_indices].get_positions()
            tm_positions = layer_atoms[tm_indices].get_positions()
            distances = get_distances(ir_positions, tm_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Find all distances within threshold
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_num}_Ir-TM'] = round(np.mean(valid_distances), 2)
        
        if len(ir_indices) > 1:
            # Calculate Ir-Ir distances within layer
            ir_positions = layer_atoms[ir_indices].get_positions()
            distances = get_distances(ir_positions, ir_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Get non-zero distances within threshold
            distances = distances[distances > 0]
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_num}_Ir-Ir'] = round(np.mean(valid_distances), 2)
        
        if len(tm_indices) > 1:
            # Calculate TM-TM distances within layer
            tm_positions = layer_atoms[tm_indices].get_positions()
            distances = get_distances(tm_positions, tm_positions, cell=atoms.get_cell(), pbc=True)[1]
            # Get non-zero distances within threshold
            distances = distances[distances > 0]
            valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
            if len(valid_distances) > 0:
                layer_distances[f'L{layer_num}_TM-TM'] = round(np.mean(valid_distances), 2)
        
        # Calculate distances between layers
        if layer_idx < len(layer_indices) - 1:
            next_layer_num = layer_num + 1
            next_layer_atoms = sorted_atoms[layer_indices[layer_idx + 1]]
            next_ir_indices = [i for i, atom in enumerate(next_layer_atoms) if atom.symbol == 'Ir']
            next_tm_indices = [i for i, atom in enumerate(next_layer_atoms) if atom.symbol in ['Mn', 'Fe', 'Co', 'Ni']]
            
            # Calculate Ir-TM distances between layers
            if ir_indices and next_tm_indices:
                ir_positions = layer_atoms[ir_indices].get_positions()
                next_tm_positions = next_layer_atoms[next_tm_indices].get_positions()
                distances = get_distances(ir_positions, next_tm_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_num}-L{next_layer_num}_Ir-TM'] = round(np.mean(valid_distances), 2)
            
            if tm_indices and next_ir_indices:
                tm_positions = layer_atoms[tm_indices].get_positions()
                next_ir_positions = next_layer_atoms[next_ir_indices].get_positions()
                distances = get_distances(tm_positions, next_ir_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_num}-L{next_layer_num}_TM-Ir'] = round(np.mean(valid_distances), 2)
            
            if ir_indices and next_ir_indices:
                # Calculate Ir-Ir distances between layers
                ir_positions = layer_atoms[ir_indices].get_positions()
                next_ir_positions = next_layer_atoms[next_ir_indices].get_positions()
                distances = get_distances(ir_positions, next_ir_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_num}-L{next_layer_num}_Ir-Ir'] = round(np.mean(valid_distances), 2)
            
            if tm_indices and next_tm_indices:
                # Calculate TM-TM distances between layers
                tm_positions = layer_atoms[tm_indices].get_positions()
                next_tm_positions = next_layer_atoms[next_tm_indices].get_positions()
                distances = get_distances(tm_positions, next_tm_positions, cell=atoms.get_cell(), pbc=True)[1]
                valid_distances = distances[distances < 3.0]  # Threshold for neighboring atoms
                if len(valid_distances) > 0:
                    layer_distances[f'L{layer_num}-L{next_layer_num}_TM-TM'] = round(np.mean(valid_distances), 2)
    
    # 인접한 레이어 간 z좌표 차이 계산
    for i in range(1, len(layer_indices)):
        current_layer = i + 1
        prev_layer = i
        if f'L{current_layer}_z' in layer_z_avg and f'L{prev_layer}_z' in layer_z_avg:
            layer_distances[f'L{prev_layer}-L{current_layer}_z_diff'] = round(
                abs(layer_z_avg[f'L{current_layer}_z'] - layer_z_avg[f'L{prev_layer}_z']), 2)
    
    return layer_charges, layer_distances, layer_z_avg

# 가장 안정한 site만 남기기
best_results = dict()  # key: (adsorbate, coverage, system_name)
for json_file in json_files:
    path_parts = json_file.split('/')
    adsorbate = path_parts[-4]
    system_name = path_parts[-3]
    site = path_parts[-2]
    # coverage 추출
    coverage = 'layer' if 'layer' in site else ('atom' if 'atom' in site else None)
    if coverage is None:
        continue
    # final_with_calculator.json 경로
    final_json = os.path.join(os.path.dirname(json_file), 'final_with_calculator.json')
    if not os.path.exists(final_json):
        continue
    try:
        atoms = read(final_json)
        energy = atoms.get_potential_energy()
    except:
        continue
    key = (adsorbate, coverage, system_name)
    if (key not in best_results) or (energy < best_results[key]['energy']):
        best_results[key] = {
            'energy': energy,
            'json_file': json_file,
            'site': site
        }

# results에 best_results만 추가
results = []
for (adsorbate, coverage, system_name), info in best_results.items():
    atoms = read(info['json_file'])
    layer_charges, layer_distances, layer_z_avg = get_layer_charges_and_distances(atoms)
    results.append({
        'Adsorbate': adsorbate,
        'System': system_name,
        'Coverage': coverage,
        'Site': info['site'],
        **layer_charges,
        **layer_distances,
        **layer_z_avg
    })

# slab 데이터는 기존 방식 유지
for json_file in slab_json_files:
    # Extract path information
    path_parts = json_file.split('/')
    system_name = path_parts[-2]  # e.g., IrFe3_111
    
    # Read the atoms with Bader charges
    atoms = read(json_file)
    
    # Get layer charges and distances based on z-coordinates
    layer_charges, layer_distances, layer_z_avg = get_layer_charges_and_distances(atoms)
    
    # Add to results
    result = {
        'Adsorbate': 'slab',
        'System': system_name,
        'Site': 'slab',
        **layer_charges,
        **layer_distances,
        **layer_z_avg
    }
    results.append(result)

# Convert to DataFrame
df = pd.DataFrame(results)

# 그래프에 사용한 데이터만 저장
df_graph = df[df['Adsorbate'] != 'slab'].copy()
graph_tsv_output = os.path.join(base_dir, 'figures', 'bader_charges_for_graph.tsv')
graph_csv_output = os.path.join(base_dir, 'figures', 'bader_charges_for_graph.csv')
df_graph.to_csv(graph_tsv_output, sep='\t', index=False, float_format='%.2f')
df_graph.to_csv(graph_csv_output, index=False, float_format='%.4f')
print(f"Graph data saved to {graph_tsv_output} and {graph_csv_output}")

# 기존 전체 데이터 저장
tsv_output = os.path.join(base_dir, 'figures', 'bader_charges_and_distances.tsv')
df.to_csv(tsv_output, sep='\t', index=False, float_format='%.2f')
csv_output = os.path.join(base_dir, 'figures', 'bader_charges_and_distances.csv')
df.to_csv(csv_output, index=False, float_format='%.4f')
print(f"Results saved to {tsv_output} and {csv_output}")

# 그래프: 흡착물+커버리지별로 하나씩만 그림
sns.set_theme(style="ticks")
for ads in selected_adsorbates:
    for cov in coverages:
        subset = df[(df['Adsorbate'] == ads) & (df['Coverage'] == cov)].copy()  # copy() 추가로 SettingWithCopyWarning 방지
        if len(subset) == 0:
            continue
        
        # 원하는 순서대로 정렬 (selected_systems의 순서 사용)
        system_order = {sys: i for i, sys in enumerate(selected_systems)}
        subset['Order'] = subset['System'].map(system_order)
        subset = subset.sort_values('Order').reset_index(drop=True)
        
        # 순서대로 정렬된 시스템 목록
        available_systems = subset['System'].tolist()
        # 시스템 이름 매핑 (0_Ir -> Ir 형식으로 변환)
        system_name_mapping = {sys: name for sys, name in zip(selected_systems, system_names)}
        x = range(len(available_systems))
        
        # Bader charge plot
        plt.figure(figsize=(4, 3))
        blues = plt.cm.Blues(np.linspace(0.9, 0.4, 4))  # Ir용 블루 컬러맵
        reds = plt.cm.Reds(np.linspace(0.9, 0.4, 4))    # TM용 레드 컬러맵
        
        # 금속, 레이어 순으로 정렬하기 위해 루프 순서 변경
        # 먼저 Ir의 모든 레이어 (L4부터 L1 순서로)
        for i, layer in zip(range(4), range(4, 0, -1)):
            if f'L{layer}_Ir' in subset.columns and not subset[f'L{layer}_Ir'].isna().all():
                charges = subset[f'L{layer}_Ir']
                plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir L{layer}')
        
        # 그 다음 TM의 모든 레이어 (L4부터 L1 순서로)
        for i, layer in zip(range(4), range(4, 0, -1)):
            if f'L{layer}_TM' in subset.columns and not subset[f'L{layer}_TM'].isna().all():
                charges = subset[f'L{layer}_TM']
                plt.plot(x, charges, marker='s', color=reds[i], label=f'TM L{layer}')
        
        plt.ylabel('Bader Charge')
        # 시스템 코드명(0_Ir 등)을 표시 이름(Ir 등)으로 변환
        display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
        plt.xticks(x, display_names)
        plt.yticks(np.arange(-1.25, 1.01, 0.25))
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plot_output = os.path.join(base_dir, 'figures', f'bader_charges_{ads}_{cov}.png')
        os.makedirs(os.path.dirname(plot_output), exist_ok=True)
        plt.savefig(plot_output, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to {plot_output}")
        
        # 레이어 간 거리 그래프는 제거 (comparison 그래프로 대체)

# slab (clean surface) 데이터에 대한 그래프 생성
slab_subset = df[df['Adsorbate'] == 'slab'].copy()
if len(slab_subset) > 0:
    # 원하는 순서대로 정렬 (selected_systems의 순서 사용)
    system_order = {sys: i for i, sys in enumerate(selected_systems)}
    slab_subset['Order'] = slab_subset['System'].map(system_order)
    slab_subset = slab_subset.sort_values('Order').reset_index(drop=True)
    
    # 순서대로 정렬된 시스템 목록
    available_systems = slab_subset['System'].tolist()
    # 시스템 이름 매핑 (0_Ir -> Ir 형식으로 변환)
    system_name_mapping = {sys: name for sys, name in zip(selected_systems, system_names)}
    x = range(len(available_systems))
    
    # Bader charge plot
    plt.figure(figsize=(4, 3))
    blues = plt.cm.Blues(np.linspace(0.9, 0.4, 4))  # Ir용 블루 컬러맵
    reds = plt.cm.Reds(np.linspace(0.9, 0.4, 4))    # TM용 레드 컬러맵
    
    # 금속, 레이어 순으로 정렬하기 위해 루프 순서 변경
    # 먼저 Ir의 모든 레이어 (L4부터 L1 순서로)
    for i, layer in zip(range(4), range(4, 0, -1)):
        if f'L{layer}_Ir' in slab_subset.columns and not slab_subset[f'L{layer}_Ir'].isna().all():
            charges = slab_subset[f'L{layer}_Ir']
            plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir L{layer}')
    
    # 그 다음 TM의 모든 레이어 (L4부터 L1 순서로)
    for i, layer in zip(range(4), range(4, 0, -1)):
        if f'L{layer}_TM' in slab_subset.columns and not slab_subset[f'L{layer}_TM'].isna().all():
            charges = slab_subset[f'L{layer}_TM']
            plt.plot(x, charges, marker='s', color=reds[i], label=f'TM L{layer}')
    
    plt.ylabel('Bader Charge')
    # 시스템 코드명(0_Ir 등)을 표시 이름(Ir 등)으로 변환
    display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
    plt.xticks(x, display_names)
    plt.yticks(np.arange(-1.25, 1.01, 0.25))
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plot_output = os.path.join(base_dir, 'figures', 'bader_charges_slab.png')
    os.makedirs(os.path.dirname(plot_output), exist_ok=True)
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}")
    
    # 레이어 간 거리 그래프는 제거 (comparison 그래프로 대체)

# ============== 새로운 그래프: L3-L4 레이어 사이 거리 - 흡착물을 범례로 ================
# 새로운 그래프 생성 - 표면 종류를 x축으로, 흡착물을 범례로 사용

# 각 레이어 간 거리에 대한 비교 그래프를 위한 함수
def create_layer_distance_comparison(layer_pair, df, selected_systems, system_names, selected_adsorbates, coverages, base_dir):
    prev_layer, current_layer = layer_pair
    col_name = f'L{prev_layer}-L{current_layer}_z_diff'
    
    # 모든 표면에 대해 하나의 그래프 생성
    plt.figure(figsize=(6, 4))
    
    # 색상 맵 정의
    colors = {
        'slab': 'black',
        '1_H': 'blue',
        '2_OH': 'green',
        '3_O': 'red'
    }
    
    # 불투명도 정의
    alpha_values = {
        'layer': 1.0,
        'atom': 0.5
    }
    
    # 전체 데이터에서 시스템 순서대로 정렬
    system_order = {sys: i for i, sys in enumerate(selected_systems)}
    system_name_mapping = {sys: name for sys, name in zip(selected_systems, system_names)}
    
    # 모든 흡착물에 대한 데이터 하나의 그래프에 표시
    for adsorbate in ['slab'] + selected_adsorbates:
        if adsorbate == 'slab':
            # slab 데이터
            subset = df[df['Adsorbate'] == 'slab'].copy()
            if len(subset) > 0:
                subset['Order'] = subset['System'].map(system_order)
                subset = subset.sort_values('Order').reset_index(drop=True)
                
                # 레이어 간 거리가 존재하는 경우에만 그래프 표시
                if col_name in subset.columns and not subset[col_name].isna().all():
                    available_systems = subset['System'].tolist()
                    display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
                    x = range(len(available_systems))
                    
                    distances = subset[col_name]
                    plt.plot(x, distances, marker='o', linestyle='-', 
                             color=colors.get(adsorbate, 'gray'), 
                             label=f'Clean Surface')
        else:
            # 흡착물 데이터 - 각 coverage별로 표시
            for cov in coverages:
                subset = df[(df['Adsorbate'] == adsorbate) & (df['Coverage'] == cov)].copy()
                if len(subset) > 0:
                    subset['Order'] = subset['System'].map(system_order)
                    subset = subset.sort_values('Order').reset_index(drop=True)
                    
                    # 레이어 간 거리가 존재하는 경우에만 그래프 표시
                    if col_name in subset.columns and not subset[col_name].isna().all():
                        available_systems = subset['System'].tolist()
                        display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
                        x = range(len(available_systems))
                        
                        # 흡착물 표시 이름 설정
                        if adsorbate == '1_H':
                            ads_name = 'H'
                        elif adsorbate == '2_OH':
                            ads_name = 'OH'
                        elif adsorbate == '3_O':
                            ads_name = 'O'
                        else:
                            ads_name = adsorbate
                        
                        distances = subset[col_name]
                        plt.plot(x, distances, marker='o', linestyle='-', 
                                 color=colors.get(adsorbate, 'gray'),
                                 alpha=alpha_values.get(cov, 1.0),
                                 label=f'{ads_name} ({cov})')
    
    # x축 설정 - 모든 가능한 시스템 표시
    all_systems = []
    for sys_code, sys_name in zip(selected_systems, system_names):
        if df['System'].str.contains(sys_code).any():
            all_systems.append((sys_code, sys_name))
    
    x_ticks = range(len(all_systems))
    x_labels = [name for _, name in all_systems]
    
    plt.xticks(x_ticks, x_labels)
    # x축 범위 설정 (-0.5부터 4.5까지)
    plt.xlim(-0.5, len(all_systems) - 0.5)
    plt.ylabel(f'L{prev_layer}-L{current_layer} Distance (Å)')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plot_output = os.path.join(base_dir, 'figures', f'layer_distances_L{prev_layer}-L{current_layer}_comparison_layer.png')
    os.makedirs(os.path.dirname(plot_output), exist_ok=True)
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}")

# L3-L4 레이어 간 거리 비교 그래프
create_layer_distance_comparison((3, 4), df, selected_systems, system_names, selected_adsorbates, coverages, base_dir)

# L2-L3 레이어 간 거리 비교 그래프
create_layer_distance_comparison((2, 3), df, selected_systems, system_names, selected_adsorbates, coverages, base_dir)

# L1-L2 레이어 간 거리 비교 그래프
create_layer_distance_comparison((1, 2), df, selected_systems, system_names, selected_adsorbates, coverages, base_dir) 