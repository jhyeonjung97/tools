from ase.io import read
import os
import pandas as pd
import matplotlib.pyplot as plt
from math import gcd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--oxygen', type=float, default=-4.5, help='Oxygen chemical potential (eV/atom)')
parser.add_argument('-x', '--fig-width', type=float, default=6, help='Figure width (inches)')
parser.add_argument('-y', '--fig-height', type=float, default=4, help='Figure height (inches)')
args = parser.parse_args()
oxygen = args.oxygen
fig_width = args.fig_width
fig_height = args.fig_height

base_path = "/Users/jiuy97/Desktop/3_RuO2/4_high_valence"
elements = ['Ru', 'Hf', 'Ta', 'W', 'Re', 'Os']

dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
dirs.sort()

# 데이터를 딕셔너리로 수집 (main_dir -> sub_dir -> atoms)
data_dict = {}

def get_element_from_sub_dir(sub_dir):
    # sub_dir 형식: '0_Ru', '1_Hf' 등
    parts = sub_dir.split('_')
    if len(parts) >= 2:
        element = '_'.join(parts[1:])  # 'RuO4', 'HfO2' 같은 경우도 처리
        # elements 리스트에 있는 것만 사용 (예: 'RuO4' -> 'Ru'로 매칭)
        for el in elements:
            if element.startswith(el):
                return el
    return None

def get_simplest_formula(element_count, oxygen_count):
    """원자 개수로부터 가장 간단한 화학식을 구하는 함수
    예: element_count=8, oxygen_count=20 -> (2, 5) -> 'M2O5'
    """
    if element_count == 0 or oxygen_count == 0:
        return None, None
    
    # 최대공약수 구하기
    common_divisor = gcd(element_count, oxygen_count)
    
    # 최대공약수로 나누어 가장 간단한 정수 비 구하기
    element_ratio = element_count // common_divisor
    oxygen_ratio = oxygen_count // common_divisor
    
    return element_ratio, oxygen_ratio

# 모든 메인 폴더 순회
for main_dir in dirs:
    main_path = os.path.join(base_path, main_dir)
    
    # 각 메인 폴더의 세부 폴더들 가져오기
    sub_dirs = [d for d in os.listdir(main_path) if os.path.isdir(os.path.join(main_path, d))]
    sub_dirs.sort()
    
    # main_dir을 키로 하는 딕셔너리 초기화
    data_dict[main_dir] = {}
    
    # 각 세부 폴더 순회
    for sub_dir in sub_dirs:
        sub_path = os.path.join(main_path, sub_dir)
        element = get_element_from_sub_dir(sub_dir)
        if element is None:
            continue
        file_path = os.path.join(sub_path, 'final_with_calculator.json')
        if os.path.exists(file_path):
            atoms = read(file_path)
            data_dict[main_dir][element] = atoms
        else:
            data_dict[main_dir][element] = None

# DataFrame 생성
df_data = {}
for main_dir in dirs:
    df_data[main_dir] = [data_dict.get(main_dir, {}).get(element, None) for element in elements]

df = pd.DataFrame(df_data, index=elements)

plt.figure(figsize=(fig_width, fig_height))
atoms_ruo2 = data_dict[dirs[0]]['Ru']
energy_ruo2 = atoms_ruo2.get_potential_energy()
plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green')
plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green')
plt.text(0.0, 0.005, 'RuO$_2$', ha='center', va='bottom')
plt.text(1.0, 0.005, 'MO$_2$', ha='center', va='bottom')
plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=0)
cmap = plt.colormaps['YlOrRd']
non_ru_elements = [el for el in elements if el != 'Ru']
colors = [cmap(i / (len(non_ru_elements) - 1)) for i in range(len(non_ru_elements))]
for i, element in enumerate(non_ru_elements):
    atoms_mo2 = data_dict[dirs[1]][element]
    atoms_mruo2 = data_dict[dirs[0]][element]
    formation_energy = (
        8*atoms_mruo2.get_potential_energy() 
        - 7*energy_ruo2 
        - 1*atoms_mo2.get_potential_energy()
    )/8/8
    plt.scatter(7/8, formation_energy, marker='D', edgecolor='black', facecolor=colors[i])
    plt.text(7/8-0.02, formation_energy, f'{element}-RuO$_2$', ha='right', va='center')
plt.xlabel(r'x, Ru$_x$M$_{1-x}$O$_2$')
plt.ylabel('Formation energy (eV/unit)')
plt.xlim(-0.1, 1.1)
# plt.ylim(-0.16, 0.14)
plt.ylim(-0.1, 0.4)
plt.savefig('RuO2_MO2_convex_hull.png', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
plt.close()

for element in non_ru_elements:
    plt.figure(figsize=(fig_width, fig_height))
    plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green')
    plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green')
    plt.text(0.0, 0.02, 'RuO$_2$', ha='center', va='bottom')
    plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=0)
    atoms_mruo2 = data_dict[dirs[0]][element]
    atoms_mxoy = data_dict[dirs[2]][element]
    element_count = atoms_mxoy.get_chemical_symbols().count(element)
    oxygen_count = atoms_mxoy.get_chemical_symbols().count('O')
    element_ratio, oxygen_ratio = get_simplest_formula(element_count, oxygen_count)
    if element_ratio == 1:
        formula = f'{element}O$_{oxygen_ratio}$'
    else:
        formula = f'{element}$_{element_ratio}$O$_{oxygen_ratio}$'
    plt.text(1.0, 0.02, formula, ha='center', va='bottom')
    
    formation_energy = (
        atoms_mruo2.get_potential_energy() 
        - energy_ruo2*7/8
        - atoms_mxoy.get_potential_energy()/element_count
        - oxygen*(2-oxygen_count/element_count)
    )/8
    plt.scatter(7/8, formation_energy, marker='D', edgecolor='black', facecolor='orange')
    plt.text(7/8-0.02, formation_energy, f'{element}-RuO$_2$', ha='right', va='center')
    plt.xlabel(r'x, Ru$_x$M$_{1-x}$O$_2$')
    plt.ylabel('Formation energy (eV/unit)')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 0.4)
    plt.savefig(f'{element}-RuO2_MO2_convex_hull.png', dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()
    plt.close()