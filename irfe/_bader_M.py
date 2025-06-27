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
selected_systems = ['1_Mn', '2_Fe', '3_Co', '4_Ni', '0_Ir']
system_names = ['Mn', 'Fe', 'Co', 'Ni', 'Ir']

# 사용할 흡착물 리스트
selected_adsorbates = ['1_H', '2_OH', '3_O']

# 커버리지 설정 (pourbaix.py와 동일)
coverage_settings = {
    '1_H': '2_layer_hol',
    '2_OH': '2_layer_brg',
    '3_O': '2_layer_hol'
}

# 0_Ir에 대해서만 1_H의 커버리지를 1_layer_top으로 변경
def get_coverage(ads, system):
    if system == '0_Ir' and ads == '1_H':
        return '1_layer_top'
    return coverage_settings[ads]

# Find all atoms_bader_charge.json files
pattern = os.path.join(base_dir, '*_*/*_*/*_*_*/atoms_bader_charge.json')
json_files = [f for f in sorted(glob.glob(pattern))
              if any(f'/{sys}/' in f for sys in selected_systems)
              and any(f'/{ads}/' in f for ads in selected_adsorbates)]

# Find slab atoms_bader_charge.json files
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
        layer_num = layer_idx + 1
        layer_atoms = sorted_atoms[indices]
        layer_z_positions = layer_atoms.get_positions()[:, 2]
        layer_z_avg[f'L{layer_num}_z'] = round(np.mean(layer_z_positions), 2)
        
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
    
    # 인접한 레이어 간 z좌표 차이 계산
    for i in range(1, len(layer_indices)):
        current_layer = i + 1
        prev_layer = i
        if f'L{current_layer}_z' in layer_z_avg and f'L{prev_layer}_z' in layer_z_avg:
            layer_distances[f'L{prev_layer}-L{current_layer}'] = round(
                abs(layer_z_avg[f'L{current_layer}_z'] - layer_z_avg[f'L{prev_layer}_z']), 2)
    
    return layer_charges, layer_distances, layer_z_avg

# 가장 안정한 site만 남기기
best_results = dict()  # key: (adsorbate, system_name)
for json_file in json_files:
    path_parts = json_file.split('/')
    adsorbate = path_parts[-4]
    system_name = path_parts[-3]
    site = path_parts[-2]
    
    # 원하는 커버리지인지 확인
    expected_coverage = get_coverage(adsorbate, system_name)
    if expected_coverage not in site:
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
        'Adsorbate': adsorbate,
        'System': system_name,
        'Site': info['site'],
        **layer_charges,
        **layer_distances,
        **layer_z_avg
    })

# slab 데이터 추가
for json_file in slab_json_files:
    path_parts = json_file.split('/')
    system_name = path_parts[-2]
    
    try:
        atoms = read(json_file)
        layer_charges, layer_distances, layer_z_avg = get_layer_charges_and_distances(atoms)
        
        result = {
            'Adsorbate': 'slab',
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

# 데이터 저장
tsv_output = os.path.join(base_dir, 'figures', 'bader_charges.tsv')
csv_output = os.path.join(base_dir, 'figures', 'bader_charges.csv')
df.to_csv(tsv_output, sep='\t', index=False, float_format='%.2f')
df.to_csv(csv_output, index=False, float_format='%.4f')
print(f"Results saved to {tsv_output} and {csv_output}")

# 그래프 생성
sns.set_theme(style="ticks")
for ads in ['slab'] + selected_adsorbates:
    subset = df[df['Adsorbate'] == ads].copy()
    if len(subset) == 0:
        continue
    
    # 원하는 순서대로 정렬
    system_order = {sys: i for i, sys in enumerate(selected_systems)}
    subset['Order'] = subset['System'].map(system_order)
    subset = subset.sort_values('Order').reset_index(drop=True)
    
    # 순서대로 정렬된 시스템 목록
    available_systems = subset['System'].tolist()
    system_name_mapping = {sys: name for sys, name in zip(selected_systems, system_names)}
    x = range(len(available_systems))
    
    # Bader charge plot
    plt.figure(figsize=(4, 3))
    blues = plt.cm.Blues(np.linspace(0.9, 0.4, 4))  # Ir용 블루 컬러맵
    reds = plt.cm.Reds(np.linspace(0.9, 0.4, 4))    # TM용 레드 컬러맵
    
    # Ir의 모든 레이어 (L4부터 L1 순서로)
    for i, layer in zip(range(4), range(4, 0, -1)):
        if f'L{layer}_Ir' in subset.columns and not subset[f'L{layer}_Ir'].isna().all():
            charges = subset[f'L{layer}_Ir']
            plt.plot(x, charges, marker='o', color=blues[i], label=f'Ir L{layer}')
    
    # TM의 모든 레이어 (L4부터 L1 순서로)
    for i, layer in zip(range(4), range(4, 0, -1)):
        if f'L{layer}_TM' in subset.columns and not subset[f'L{layer}_TM'].isna().all():
            charges = subset[f'L{layer}_TM']
            plt.plot(x, charges, marker='s', color=reds[i], label=f'TM L{layer}')
    
    plt.ylabel('Bader Charge (e$^{\mathrm{-}}$)')
    display_names = [system_name_mapping.get(sys, sys) for sys in available_systems]
    plt.xticks(x, display_names)
    plt.yticks(np.arange(-1.25, 1.01, 0.25))
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', handletextpad=0.5, labelspacing=0.2)
    
    # 흡착물 이름 설정
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
    
    plot_output = os.path.join(base_dir, 'figures', f'bader_charges_{ads_name}.png')
    os.makedirs(os.path.dirname(plot_output), exist_ok=True)
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}")

# 레이어 간 거리 그래프
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
    slab_subset = df[df['Adsorbate'] == 'slab'].copy()
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
        subset = df[df['Adsorbate'] == ads].copy()
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
    
    plot_output = os.path.join(base_dir, 'figures', f'layer_distances_L{prev_layer}-L{current_layer}.png')
    plt.savefig(plot_output, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_output}") 