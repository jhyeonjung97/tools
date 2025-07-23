import os
import re
import socket
import argparse
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
import glob
from ase.io import read
from mendeleev import element
from collections import Counter
import pandas as pd
import os
# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--gc', action='store_true', help='Enable GC-DFT mode')
parser.add_argument('--bulk', action='store_true', help='Enable bulk Pourbaix mode')
parser.add_argument('--suffix', type=str, default='', help='Suffix for output filename')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--save-dir', action='store_true', help='Save to predefined directory')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--figx', type=float, default=4, help='Figure width in inches (default: 6)')
parser.add_argument('--figy', type=float, default=3, help='Figure height in inches (default: 7)')
parser.add_argument('--json-dir', type=str, default='.', help='json 파일이 있는 폴더 경로 (기본값: 현재 폴더)')
parser.add_argument('--csv-dir', type=str, default='.', help='label.csv가 있는 폴더 경로 (기본값: 현재 폴더)')
args = parser.parse_args()

def main():
    elements = set()
    json_dir = args.json_dir
    json_files = glob.glob(os.path.join(json_dir, "*.json"))

    csv_dir = args.csv_dir
    label_csv_path = os.path.join(csv_dir, 'label.csv')
    
    file_labels = {}
    if os.path.exists(label_csv_path):
        label_df = pd.read_csv(label_csv_path, header=None)
        for idx, row in label_df.iterrows():
            file_labels[row[0]] = row[1]
    else:
        for json_file in json_files:
            atoms = read(json_file)
            formula = atoms.get_chemical_formula()
            file_labels[os.path.basename(json_file)] = formula
    # 원자번호 순으로 정렬
    for json_file in json_files:
        atoms = read(json_file)
        elements.update(atoms.get_chemical_symbols())
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    print(f"\n전체 구성 원소(원자번호 순): {', '.join(sorted_elements)}")
    print("\n파일별 라벨:")
    for fname, label in file_labels.items():
        print(f"{fname}: {label}")

    # surfs DataFrame 생성
    # 원자번호 순 전체 원소 리스트
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    # 각 원소별 최소 개수 저장용 딕셔너리 (초기값은 None)
    min_counts = {el: None for el in sorted_elements}
    surfs = []
    for json_file in json_files:
        atoms = read(json_file)
        energy = atoms.get_potential_energy()
        symbols = atoms.get_chemical_symbols()
        row = {}
        row['E_DFT'] = energy
        row['e'] = 0
        # 각 원소별 개수 세기
        symbol_count = Counter(symbols)
        for el in sorted_elements:
            count = symbol_count[el]
            row[el] = count
            # 최소값 갱신
            if min_counts[el] is None:
                min_counts[el] = count
            else:
                min_counts[el] = min(min_counts[el], count)
        # 파일명 또는 라벨을 name에 저장
        row['name'] = file_labels.get(os.path.basename(json_file), json_file)
        surfs.append(row)

    print("\n각 원소별 모든 파일에서의 최소 개수:")
    for el in sorted_elements:
        print(f"{el}: {min_counts[el]}")

    # 각 row에서 최소 개수만큼 빼기
    for row in surfs:
        for el in sorted_elements:
            row[el] = row[el] - min_counts[el]
    
    # 모든 값이 0인 원소 컬럼 제거
    zero_cols = [el for el in sorted_elements if all(row[el] == 0 for row in surfs)]
    remaining_elements = [el for el in sorted_elements if el not in zero_cols]
    
    # surfs에서 zero_cols에 해당하는 키 삭제
    for row in surfs:
        for col in zero_cols:
            del row[col]
    
    # O, H를 제외한 원소들을 unique_elements로 정의
    unique_elements = [el for el in remaining_elements if el not in ['O', 'H']]
    print("\nsorted_elements:", sorted_elements)
    print("remaining_elements:", remaining_elements)
    print("unique_elements:", unique_elements)
    
    print("\n", surfs)
    surfs_df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['name'])
    print("\n", surfs_df)


if __name__ == "__main__":
    main()