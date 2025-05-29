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

# 사용할 흡착물 리스트 - layer_top만 처리
selected_adsorbates = ['1_H', '2_OH', '3_O']

# Find all atoms_bader_charge.json files for layer_top data
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = [f for f in sorted(glob.glob(pattern))
              if any(f'/{sys}/' in f for sys in selected_systems)
              and any(f'/{ads}/' in f for ads in selected_adsorbates)
              and 'layer_top' in f]

# Find slab atoms_bader_charge.json files (선태된 표면만)
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
    
    # 모든 금속 원자 간 거리 계산을 위한 리스트
    all_metal_distances = []
    
    # 각 레이어별 결합 수를 저장할 딕셔너리
    layer_bonds = {}
    
    # 모든 금속 원자의 인덱스와 위치를 저장할 리스트
    all_metal_indices = []
    all_metal_positions = []
    
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
        
        # 현재 레이어의 금속 원자 인덱스와 위치 저장
        metal_indices = ir_indices + tm_indices
        if metal_indices:
            all_metal_indices.extend([indices[i] for i in metal_indices])
            all_metal_positions.extend(layer_atoms[metal_indices].get_positions())
    
    # 전체 셀의 모든 금속 원자 간 거리 계산
    if len(all_metal_positions) > 1:
        print("\n전체 금속 원자 수:", len(all_metal_indices))
        all_metal_positions = np.array(all_metal_positions)
        distances = get_distances(all_metal_positions, all_metal_positions, cell=atoms.get_cell(), pbc=True)[1]
        # 0이 아닌 거리만 선택 (같은 원자 간 거리 제외)
        valid_distances = distances[distances > 0]
        print("총 결합 수:", len(valid_distances))
        # 3Å 이내의 거리만 선택
        # valid_distances = valid_distances[valid_distances < 3.0]
        print("고유한 결합 수 (중복 제거):", len(valid_distances))
        all_metal_distances.extend(valid_distances)
    
    # 인접한 레이어 간 z좌표 차이 계산
    for i in range(1, len(layer_indices)):
        current_layer = i + 1
        prev_layer = i
        if f'L{current_layer}_z' in layer_z_avg and f'L{prev_layer}_z' in layer_z_avg:
            layer_distances[f'L{prev_layer}-L{current_layer}'] = round(
                abs(layer_z_avg[f'L{current_layer}_z'] - layer_z_avg[f'L{prev_layer}_z']), 2)
    
    # 모든 금속 원자 간 거리의 평균 계산
    if all_metal_distances:
        layer_distances['avg'] = round(np.mean(all_metal_distances), 2)
    
    print(np.mean(all_metal_distances))
    
    return layer_charges, layer_distances, layer_z_avg

# 가장 안정한 layer_top site만 남기기
best_results = dict()  # key: (adsorbate, system_name)
for json_file in json_files:
    if 'layer_top' not in json_file:
        continue
        
    path_parts = json_file.split('/')
    adsorbate = path_parts[-4]
    system_name = path_parts[-3]
    site = path_parts[-2]
    
    # final_with_calculator.json 경로
    final_json = os.path.join(os.path.dirname(json_file), 'final_with_calculator.json')
    if not os.path.exists(final_json):
        continue
    
    try:
        atoms = read(final_json)
        energy = atoms.get_potential_energy()
    except:
        continue
    
    key = (adsorbate, system_name)
    if (key not in best_results) or (energy < best_results[key]['energy']):
        best_results[key] = {
            'energy': energy,
            'json_file': json_file,
            'site': site
        }

# results에 best_results만 추가
for (adsorbate, system_name), info in best_results.items():
    atoms = read(info['json_file'])
    layer_charges, layer_distances, layer_z_avg = get_layer_charges_and_distances(atoms)
    results.append({
        'Ads': adsorbate,
        'System': system_name,
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
    try:
        atoms = read(json_file)
        
        # Get layer charges and distances based on z-coordinates
        layer_charges, layer_distances, layer_z_avg = get_layer_charges_and_distances(atoms)
        
        # Add to results
        result = {
            'Ads': 'slab',
            'System': system_name,
            'Site': 'slab',
            **layer_charges,
            **layer_distances,
            **layer_z_avg
        }
        results.append(result)
    except:
        print(f"Warning: Could not process {json_file}")
        continue

# Convert to DataFrame
df = pd.DataFrame(results)

# Adsorbate를 Ads로 변경
df = df.rename(columns={'Ads': 'Ads'})

# 필요한 컬럼만 선택
needed_columns = [
    'Ads', 'System',
    'L1_Ir', 'L1_TM', 'L2_Ir', 'L2_TM', 'L3_Ir', 'L3_TM', 'L4_Ir', 'L4_TM',
    'L1-L2', 'L2-L3', 'L3-L4',
    'avg',
    'L1_z', 'L2_z', 'L3_z', 'L4_z'
]
df = df[needed_columns]

# 그래프에 사용한 데이터 저장
graph_tsv_output = os.path.join(base_dir, 'bader_charges_final.tsv')
graph_csv_output = os.path.join(base_dir, 'bader_charges_final.csv')

# 금속 간 평균 거리 데이터를 포함하여 저장
df.to_csv(graph_tsv_output, sep='\t', index=False, float_format='%.2f')
df.to_csv(graph_csv_output, index=False, float_format='%.4f')
print(f"Graph data saved to {graph_tsv_output} and {graph_csv_output}")

# 데이터 확인을 위한 출력
print("\nDataFrame columns:")
print(df.columns.tolist())
print("\nSample data:")
print(df[['Ads', 'System', 'avg']].head())

# 그래프: 흡착물별로 하나씩만 그림
sns.set_theme(style="ticks")
for ads in ['slab'] + selected_adsorbates:
    subset = df[df['Ads'] == ads].copy()
    if len(subset) == 0:
        print(f"Warning: No data found for {ads}")
        continue
    
    # 원하는 순서대로 정렬
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
    
    plt.ylabel('Bader Charge (e$^{\mathrm{-}}$)')
    # 시스템 코드명(0_Ir 등)을 표시 이름(Ir 등)으로 변환
    display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
    plt.xticks(x, display_names)
    plt.yticks(np.arange(-1.25, 1.01, 0.25))
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5, labelspacing=0.2)
    
    # 원본 흡착물 이름 추출
    if ads == 'slab':
        ads_name = 'Clean'
    elif ads == '1_H':
        ads_name = 'H'
    elif ads == '2_OH':
        ads_name = 'OH'
    elif ads == '3_O':
        ads_name = 'O'
    else:
        ads_name = ads
    
    plot_output = os.path.join(base_dir, f'bader_charges_final_{ads_name}.png')
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}")

# 레이어 간 거리 그래프 - 레이어 쌍별로 생성
for prev_layer, current_layer in [(1, 2), (2, 3), (3, 4)]:
    plt.figure(figsize=(4, 3))
    
    # 색상 정의
    colors = {
        'slab': 'black',
        '1_H': 'blue',
        '2_OH': 'green',
        '3_O': 'red'
    }
    
    # Clean surface부터 시작
    slab_subset = df[df['Ads'] == 'slab'].copy()
    if len(slab_subset) > 0:
        system_order = {sys: i for i, sys in enumerate(selected_systems)}
        slab_subset['Order'] = slab_subset['System'].map(system_order)
        slab_subset = slab_subset.sort_values('Order').reset_index(drop=True)
        
        col_name = f'L{prev_layer}-L{current_layer}'
        if col_name in slab_subset.columns and not slab_subset[col_name].isna().all():
            x = range(len(slab_subset))
            distances = slab_subset[col_name]
            plt.plot(x, distances, marker='o', color=colors['slab'], label='Clean')
    
    # 각 흡착물에 대한 데이터 플롯
    for ads in selected_adsorbates:
        subset = df[df['Ads'] == ads].copy()
        if len(subset) == 0:
            continue
        
        system_order = {sys: i for i, sys in enumerate(selected_systems)}
        subset['Order'] = subset['System'].map(system_order)
        subset = subset.sort_values('Order').reset_index(drop=True)
        
        col_name = f'L{prev_layer}-L{current_layer}'
        if col_name in subset.columns and not subset[col_name].isna().all():
            x = range(len(subset))
            distances = subset[col_name]
            
            # 흡착물 이름 설정
            if ads == '1_H':
                ads_name = 'H'
            elif ads == '2_OH':
                ads_name = 'OH'
            elif ads == '3_O':
                ads_name = 'O'
            else:
                ads_name = ads
            
            plt.plot(x, distances, marker='o', color=colors[ads], label=ads_name)
    
    plt.ylabel(f'L{prev_layer}-L{current_layer} Distance (Å)')
    plt.ylim(2.00, 2.35)
    display_names = [system_name_mapping.get(sys, sys) for sys in selected_systems]
    plt.xticks(range(len(selected_systems)), display_names)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5, labelspacing=0.2)
    
    plot_output = os.path.join(base_dir, f'layer_distances_final_L{prev_layer}-L{current_layer}.png')
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}")

# 금속 간 평균 거리 그래프
plt.figure(figsize=(4, 3))

# 색상 정의
colors = {
    'slab': 'black',
    '1_H': 'blue',
    '2_OH': 'green',
    '3_O': 'red'
}

# Clean surface부터 시작
slab_subset = df[df['Ads'] == 'slab'].copy()
if len(slab_subset) > 0:
    system_order = {sys: i for i, sys in enumerate(selected_systems)}
    slab_subset['Order'] = slab_subset['System'].map(system_order)
    slab_subset = slab_subset.sort_values('Order').reset_index(drop=True)
    
    if 'avg' in slab_subset.columns and not slab_subset['avg'].isna().all():
        x = range(len(slab_subset))
        distances = slab_subset['avg']
        plt.plot(x, distances, marker='o', color=colors['slab'], label='Clean Surface')

# 각 흡착물에 대한 데이터 플롯
for ads in selected_adsorbates:
    subset = df[df['Ads'] == ads].copy()
    if len(subset) == 0:
        continue
    
    system_order = {sys: i for i, sys in enumerate(selected_systems)}
    subset['Order'] = subset['System'].map(system_order)
    subset = subset.sort_values('Order').reset_index(drop=True)
    
    if 'avg' in subset.columns and not subset['avg'].isna().all():
        x = range(len(subset))
        distances = subset['avg']
        
        # 흡착물 이름 설정
        if ads == '1_H':
            ads_name = 'H'
        elif ads == '2_OH':
            ads_name = 'OH'
        elif ads == '3_O':
            ads_name = 'O'
        else:
            ads_name = ads
        
        plt.plot(x, distances, marker='o', color=colors[ads], label=ads_name)

plt.ylabel('Average Metal-Metal Distance (Å)')
display_names = [system_name_mapping.get(sys, sys) for sys in selected_systems]
plt.xticks(range(len(selected_systems)), display_names)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5, labelspacing=0.2)

plot_output = os.path.join(base_dir, 'layer_distances_final_avg.png')
plt.savefig(plot_output, bbox_inches='tight')
plt.close()
print(f"Plot saved to {plot_output}") 

print(df)