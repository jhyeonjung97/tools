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
from pymatgen.core.ion import Ion
import json

# 전역 상수 정의
SUBSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}
SUPERSCRIPT_NUMS = {'2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--json-dir', type=str, default='.', help='json 파일이 있는 폴더 경로 (기본값: 현재 폴더)')
parser.add_argument('--csv-dir', type=str, default='.', help='label.csv가 있는 폴더 경로 (기본값: 현재 폴더)')
parser.add_argument('--conc', type=float, default=10**-6, help='concentration (기본값: 10^-6)')
parser.add_argument('--bulk', action='store_true', help='bulk 모드')
parser.add_argument('--ads', type=str, nargs='*', default=[], help='adsorbate 원소들 (예: --ads N O)')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--pHmin', type=float, default=0, help='Minimum pH value (default: 0)')
parser.add_argument('--pHmax', type=float, default=14, help='Maximum pH value (default: 14)')
parser.add_argument('--Umin', type=float, default=-2.0, help='Minimum potential value (default: -2.0)')
parser.add_argument('--Umax', type=float, default=4.0, help='Maximum potential value (default: 4.0)')
parser.add_argument('--figx', type=float, default=4, help='Figure width in inches (default: 6)')
parser.add_argument('--figy', type=float, default=3, help='Figure height in inches (default: 7)')
parser.add_argument('--show', action='store_true', help='Show the plot')
parser.add_argument('--gibbs', action='store_true', help='Apply Gibbs free energy correction from G_corr column')
parser.add_argument('--gc', action='store_true', help='Apply Grand Canonical DFT using A, B, C columns')
args = parser.parse_args()

is_bulk = args.bulk
ads_elements = args.ads
is_gibbs = args.gibbs
is_gc = args.gc

def main():
    elements = set()
    json_dir = args.json_dir
    json_files = glob.glob(os.path.join(json_dir, "*.json"))

    csv_dir = args.csv_dir
    label_csv_path = os.path.join(csv_dir, 'label.csv')
    
    file_labels = {}
    file_oh_counts = {}
    file_gibbs_corrections = {}
    file_gc_params = {}  # A, B, C 파라미터 저장
    
    if os.path.exists(label_csv_path):
        # 헤더가 있는 CSV 파일 읽기
        label_df = pd.read_csv(label_csv_path, header=0)
        for idx, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']
            if '#OH' in row:
                file_oh_counts[json_name] = float(row['#OH'])
            if is_gibbs and 'G_corr' in row:
                file_gibbs_corrections[json_name] = float(row['G_corr'])
            if is_gc and all(col in row for col in ['A', 'B', 'C']):
                file_gc_params[json_name] = {
                    'A': float(row['A']),
                    'B': float(row['B']),
                    'C': float(row['C'])
                }
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
    # print(f"\n전체 구성 원소(원자번호 순): {', '.join(sorted_elements)}")
    # print("\n파일별 라벨:")
    # for fname, label in file_labels.items():
    #     print(f"{fname}: {label}")

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
        json_basename = os.path.basename(json_file)
        
        row = {}
        row['E_DFT'] = energy
        row['e'] = 0
        # 각 원소별 개수 세기
        symbol_count = Counter(symbols)
        for el in sorted_elements:
            count = symbol_count[el]
            row[el] = float(count)
            # 최소값 갱신
            if min_counts[el] is None:
                min_counts[el] = count
            else:
                min_counts[el] = min(min_counts[el], count)
        # 파일명 또는 라벨을 name에 저장
        row['conc'] = 1
        row['name'] = file_labels.get(json_basename, json_file)
        
        # Gibbs correction 추가
        if is_gibbs and json_basename in file_gibbs_corrections:
            row['gibbs_corr'] = file_gibbs_corrections[json_basename]
        else:
            row['gibbs_corr'] = 0.0
            
        # GC 파라미터 추가
        if is_gc and json_basename in file_gc_params:
            row['A'] = file_gc_params[json_basename]['A']
            row['B'] = file_gc_params[json_basename]['B']
            row['C'] = file_gc_params[json_basename]['C']
        else:
            row['A'] = 0.0
            row['B'] = 0.0
            row['C'] = 0.0
            
        surfs.append(row)

    # print("\n각 원소별 모든 파일에서의 최소 개수:")
    # for el in sorted_elements:
    #     print(f"{el}: {min_counts[el]}")

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

    # surfs 중에서 unique_elements가 모두 0인 것들 중 에너지가 가장 낮은 것을 ref_surf로 설정
    ref_candidates = []
    for i, row in enumerate(surfs):
        # unique_elements에 해당하는 원소가 모두 0인지 확인
        if all(row[elem] == 0.0 for elem in unique_elements if elem in row):
            ref_candidates.append((i, float(row['E_DFT'])))
    
    if ref_candidates:
        # 에너지가 가장 높은 것을 ref_surf로 선택
        ref_surf_idx, ref_energy = max(ref_candidates, key=lambda x: x[1])
        ref_surf = surfs[ref_surf_idx]
        print(f"\nReference surface (index {ref_surf_idx}): {ref_surf['name']}, E_DFT: {ref_energy:.2f} eV")
        print(f"Composition: " + ", ".join([f"{elem}: {ref_surf[elem]}" for elem in remaining_elements if elem in ref_surf]))
    else:
        print("No reference surface found (no surface with all unique_elements = 0)")

    # formation energy correction 적용
    reference_surface_energy = ref_surf['E_DFT']
    for k in range(len(surfs)):
        # OH count 가져오기
        oh_count = 0
        json_basename = None
        for fname, label in file_labels.items():
            if label == surfs[k]['name']:
                json_basename = fname
                break
        
        if json_basename and json_basename in file_oh_counts:
            oh_count = file_oh_counts[json_basename]
        
        formation_energy_correction = (
            - (surfs[k]['H'] - oh_count) * (gh - dgh)
            - (surfs[k]['O'] - oh_count) * (go - dgo)
            - oh_count * (goh - dgoh)
        )
        
        # Gibbs free energy correction 적용 (159번째 라인 근처)
        gibbs_correction = surfs[k]['gibbs_corr'] if is_gibbs else 0.0
        
        surfs[k]['E_DFT'] = surfs[k]['E_DFT'] - reference_surface_energy + formation_energy_correction + gibbs_correction
    
    # thermodynamic_data.json에서 unique_elements에 해당하는 데이터 가져오기
    thermo_data_path = os.path.join(os.path.dirname(__file__), 'thermodynamic_data.json')
    thermo_data = {}
    if os.path.exists(thermo_data_path):
        with open(thermo_data_path, 'r') as f:
            thermo_data = json.load(f)
    
    # reference_energies.json 읽기
    ref_energies_path = os.path.join(os.path.dirname(__file__), 'reference_energies.json')
    ref_energies = {}
    if os.path.exists(ref_energies_path):
        with open(ref_energies_path, 'r') as f:
            ref_energies = json.load(f)
    
    # solids, ions, gases 리스트 생성
    ions = []
    solids = []
    gases = []
    
    # unique_elements에 대한 thermodynamic 데이터 처리
    for el in unique_elements:
        if el in thermo_data:
            print(f"\n{el} thermodynamic data:")
                
            # ions 처리
            if 'ions' in thermo_data[el] and thermo_data[el]['ions'] != {}:
                print(f"  Ions reduced dict:")
                for ion_formula, energy in thermo_data[el]['ions'].items():
                    try:
                        reduced_dict = Ion.from_formula(ion_formula).to_reduced_dict
                        print(f"    {ion_formula}: {reduced_dict}, energy: {energy}")
                        
                        # energy 보정 계산
                        energy_correction = 0
                        # O 개수만큼 water 에너지 더하기
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # reference_energies에 있는 원소들의 에너지 더하기
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # row 생성 (보정된 energy 사용)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = args.conc                        
                        row['name'] = format_name(ion_formula) + '(aq)'
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # 분수를 윗첨자/아랫첨자로 표시
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']
                        
                        ions.append(row)
                    except:
                        print(f"    {ion_formula}: parsing failed, energy: {energy}")

            # solids 처리
            if 'solids' in thermo_data[el] and thermo_data[el]['solids'] != {}:
                print(f"  Solids reduced dict:")
                for solid_formula, energy in thermo_data[el]['solids'].items():
                    try:
                        reduced_dict = Ion.from_formula(solid_formula).to_reduced_dict
                        print(f"    {solid_formula}: {reduced_dict}, energy: {energy}")
                        
                        # energy 보정 계산
                        energy_correction = 0
                        # O 개수만큼 water 에너지 더하기
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # reference_energies에 있는 원소들의 에너지 더하기
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # row 생성 (보정된 energy 사용)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = 1
                        row['name'] = format_name(solid_formula) + '(s)'
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # 분수를 윗첨자/아랫첨자로 표시
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']
                        
                        solids.append(row)
                    except:
                        print(f"    {solid_formula}: parsing failed, energy: {energy}")

            # gases 처리
            if 'gases' in thermo_data[el] and thermo_data[el]['gases'] != {}:
                print(f"  Gases reduced dict:")
                for gas_formula, energy in thermo_data[el]['gases'].items():
                    try:
                        reduced_dict = Ion.from_formula(gas_formula).to_reduced_dict
                        print(f"    {gas_formula}: {reduced_dict}, energy: {energy}")
                        
                        # energy 보정 계산
                        energy_correction = 0
                        # O 개수만큼 water 에너지 더하기
                        if 'O' in reduced_dict:
                            energy_correction += water * reduced_dict['O']
                        # reference_energies에 있는 원소들의 에너지 더하기
                        for elem in remaining_elements:
                            if elem in ref_energies and elem in reduced_dict:
                                energy_correction += ref_energies[elem] * reduced_dict[elem]
                        
                        # row 생성 (보정된 energy 사용)
                        row = {'E_DFT': energy/calmol + energy_correction, 'e': 0}
                        for elem in remaining_elements:
                            row[elem] = int(reduced_dict.get(elem, 0))
                        if 'charge' in reduced_dict:
                            row['e'] = int(reduced_dict['charge'])
                        row['conc'] = args.conc
                        row['name'] = format_name(gas_formula) + '(g)'
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            # 분수를 윗첨자/아랫첨자로 표시
                            subscript_count = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
                            row['name'] = f'¹⁄{subscript_count}' + row['name']
                        
                        gases.append(row)
                    except:
                        print(f"    {gas_formula}: parsing failed, energy: {energy}")
        else:
            print(f"\n{el}: No thermodynamic data found")

    nsurfs, nions, nsolids, ngases = len(surfs), len(ions), len(solids), len(gases)
    
    # DataFrame 컬럼 정의
    base_columns = ['E_DFT', 'e'] + remaining_elements + ['conc', 'name']
    surf_columns = base_columns + ['gibbs_corr', 'A', 'B', 'C'] if (is_gibbs or is_gc) else base_columns
    
    surfs_df = pd.DataFrame(surfs, columns=surf_columns)
    ions_df = pd.DataFrame(ions, columns=base_columns)
    solids_df = pd.DataFrame(solids, columns=base_columns)
    gases_df = pd.DataFrame(gases, columns=base_columns)
    
    print()
    print(format_df_for_display(surfs_df))
    print(f"Surfs: {nsurfs} entries\n")
    print(format_df_for_display(ions_df))
    print(f"Ions: {nions} entries\n")
    print(format_df_for_display(solids_df))
    print(f"Solids: {nsolids} entries\n")
    print(format_df_for_display(gases_df))
    print(f"Gases: {ngases} entries\n")

    new_surfs = []
    
    if is_bulk:
        # 각 unique_element에 대해 처리
        for unique_elem in unique_elements:
            # surfs 중에서 해당 unique_element가 0인 것들 찾기
            for k in range(nsurfs):
                if surfs[k][unique_elem] == 0:
                    # ions에서 해당 unique_element가 1개인 화합물들과 조합
                    for ion in ions:
                        if ion[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + ion[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * ion[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + ion[key]
                            new_surfs.append(new_surf)
                    
                    # solids에서 해당 unique_element가 1개인 화합물들과 조합
                    for solid in solids:
                        if solid[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + solid[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * solid[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + solid[key]
                            new_surfs.append(new_surf)
                    
                    # gases에서 해당 unique_element가 1개인 화합물들과 조합
                    for gas in gases:
                        if gas[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + gas[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * gas[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + gas[key]
                            new_surfs.append(new_surf)
    else:
        # is_bulk이 False인 경우: 각 원소의 순수한 형태만 사용 + ads 원소들은 모든 화합물 사용
        # 각 원소의 순수한 형태 화합물들 찾기 (해당 원소만 1개, 나머지는 0개)
        ref_ions = []
        ref_solids = []
        ref_gases = []
        
        for elem in unique_elements:
            # solids에서 해당 원소만 1개이고 나머지는 0개인 것들
            for solid in solids:
                if solid[elem] == 1 and all(solid[other_elem] == 0 for other_elem in remaining_elements if other_elem != elem):
                    if solid not in ref_solids:
                        ref_solids.append(solid)
            
            # gases에서 해당 원소만 1개이고 나머지는 0개인 것들
            for gas in gases:
                if gas[elem] == 1 and all(gas[other_elem] == 0 for other_elem in remaining_elements if other_elem != elem):
                    if gas not in ref_gases:
                        ref_gases.append(gas)
        
        # ads_elements에 해당하는 화합물들도 추가 (모든 조성 허용)
        for elem in ads_elements:
            if elem in unique_elements:
                # 해당 원소를 포함한 모든 화합물들 추가
                ref_ions.extend([ion for ion in ions if elem in ion and ion[elem] > 0 and ion not in ref_ions])
                ref_solids.extend([solid for solid in solids if elem in solid and solid[elem] > 0 and solid not in ref_solids])
                ref_gases.extend([gas for gas in gases if elem in gas and gas[elem] > 0 and gas not in ref_gases])
        
        # 각 unique_element에 대해 처리
        for unique_elem in unique_elements:
            # surfs 중에서 해당 unique_element가 0인 것들 찾기
            for k in range(nsurfs):
                if surfs[k][unique_elem] == 0:
                    # ref_ions에서 해당 unique_element가 1개인 화합물들과 조합
                    for ion in ref_ions:
                        if ion[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + ion[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * ion[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + ion[key]
                            new_surfs.append(new_surf)
                    
                    # ref_solids에서 해당 unique_element가 1개인 화합물들과 조합
                    for solid in ref_solids:
                        if solid[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + solid[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * solid[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + solid[key]
                            new_surfs.append(new_surf)
                    
                    # ref_gases에서 해당 unique_element가 1개인 화합물들과 조합
                    for gas in ref_gases:
                        if gas[unique_elem] == 1:
                            new_surf = {}
                            for key in surfs[k]:
                                if key == 'name':
                                    new_surf[key] = surfs[k][key] + '+' + gas[key]
                                elif key == 'conc':
                                    new_surf[key] = surfs[k][key] * gas[key]
                                elif key in ['gibbs_corr', 'A', 'B', 'C']:
                                    new_surf[key] = surfs[k][key]  # 원래 surface 값 유지
                                else:
                                    new_surf[key] = surfs[k][key] + gas[key]
                            new_surfs.append(new_surf)

    # new_surfs를 surfs에 추가
    surfs.extend(new_surfs)
    # unique_elements가 모두 0이 아닌 것들만 남기기
    surfs = [surf for surf in surfs if not all(surf[elem] == 0 for elem in unique_elements)]
    nsurfs = len(surfs)
    df = pd.DataFrame(surfs, columns=surf_columns)
    print(format_df_for_display(df))
    print(f"After adding combinations: {nsurfs} surfs entries")

    # pH and potential grid 설정
    tick = args.tick if hasattr(args, 'tick') else 0.01
    pHmin, pHmax = args.pHmin, args.pHmax
    pHrange = np.arange(pHmin, pHmax + tick, tick)
    Umin, Umax = args.Umin, args.Umax
    Urange = np.arange(Umin, Umax + 0.06 * 14, tick)
    target_pH = args.ph if hasattr(args, 'ph') else 0

    # lowest surfaces 계산
    lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

    pHindex = 0
    for pH in pHrange:
        Uindex = 0
        for U in Urange:
            values = []
            for k in range(nsurfs):
                value = dg(surfs[k], pH, U, surfs[ref_surf_idx])
                values.append(value)
            sorted_values = sorted(range(len(values)), key=lambda k: values[k])
            lowest_surfaces[Uindex][pHindex] = sorted_values[0]
            Uindex += 1
        pHindex += 1

    # 최소 좌표 찾기
    min_coords = {}
    n_rows, n_cols = lowest_surfaces.shape

    for j in range(n_cols):
        for i in range(n_rows):
            sid = int(lowest_surfaces[i, j])
            x = pHrange[j]
            y = Urange[i]
            if sid not in min_coords:
                min_coords[sid] = (x, y)
            else:
                current_x, current_y = min_coords[sid]
                if x < current_x or (x == current_x and y < current_y):
                    min_coords[sid] = (x, y)

    for sid in sorted(min_coords):
        x, y = min_coords[sid]
        name = surfs[int(sid)]['name']
        print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}, name = {name}")

    # Pourbaix diagram 그리기
    fig, ax = plt.subplots(figsize=(args.figx, args.figy), dpi=100)
    ax.axis([pHmin, pHmax, Umin, Umax])
    ax.set_xlabel('pH', labelpad=0)
    ax.set_ylabel('E (V vs. SHE)', labelpad=-6)
    ax.tick_params(right=True, direction="in")
    plt.xticks(np.arange(pHmin, pHmax + 1, 2))

    # unique surface IDs
    unique_ids = np.unique(lowest_surfaces)

    # target_pH에서의 lowest surfaces 계산
    lowest_surfaces_pH = np.zeros(len(Urange))
    second_lowest_surfaces_pH = np.zeros(len(Urange))
    for Uindex, U in enumerate(Urange):
        values = []
        for k in range(nsurfs):
            value = dg(surfs[k], pH=target_pH, U=U, ref_surf=surfs[ref_surf_idx])
            values.append(value)
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces_pH[Uindex] = sorted_values[0]
        if len(sorted_values) > 1:
            second_lowest_surfaces_pH[Uindex] = sorted_values[1]

    unique_ids_pH = np.unique(lowest_surfaces_pH.astype(int))
    unique_second_ids_pH = np.unique(second_lowest_surfaces_pH.astype(int))

    print(f"\nMost stable surfaces at pH={target_pH}:")
    for sid in unique_ids_pH:
        name = surfs[int(sid)]['name']
        print(f"Surface {sid}: {name}")

    print(f"\nSecond most stable surfaces at pH={target_pH}:")
    for sid in unique_second_ids_pH:
        name = surfs[int(sid)]['name']
        print(f"Surface {sid}: {name}")

    # 색상 매핑 및 플롯
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_ids)))
    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(len(unique_ids) + 1) - 0.5
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # ID 매핑
    id_map = {val: idx for idx, val in enumerate(unique_ids)}
    mapped_surfaces = np.vectorize(id_map.get)(lowest_surfaces)

    # 범례 생성
    for idx, surf_id in enumerate(unique_ids):
        label = surfs[int(surf_id)]['name']
        plt.plot([], [], color=colors[idx], linewidth=5, label=label)

    # pcolormesh
    pH_grid, U = np.meshgrid(pHrange, Urange)
    plt.pcolormesh(pH_grid, U, mapped_surfaces, cmap=cmap, norm=norm)

    plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., 
               fontsize='small', ncol=1, handlelength=3, edgecolor='black')

    # 물의 안정성 영역 표시
    plt.plot(pHrange, 1.23-pHrange*const, '--', lw=1, color='mediumblue')
    plt.plot(pHrange, 0-pHrange*const, '--', lw=1, color='mediumblue')

    plt.savefig(f'pourbaix_diagram.png', dpi=300, bbox_inches='tight')
    print("Pourbaix diagram saved as pourbaix_diagram.png")

    if args.show:
        plt.show()

# dg 함수 정의 (Gibbs free energy 계산)
def dg(surf, pH, U, ref_surf):
    if is_gc:
        # Grand Canonical DFT: A*U^2 + B*U + C 형태
        surface_term = (surf['A']*(U**2) + surf['B']*U + surf['C']) - (ref_surf['A']*(U**2) + ref_surf['B']*U + ref_surf['C'])
    else:
        surface_term = surf['E_DFT'] - ref_surf['E_DFT']
    
    U_coeff = surf['H'] - 2*surf['O'] - surf['e']
    pH_coeff = surf['H'] - 2*surf['O']
    dg_value = surface_term + U_coeff*U + const*pH_coeff*pH + const*log10(surf['conc'])
    return dg_value

# 이온 이름에서 + 또는 - 개수를 세어서 윗첨자로 변환하고, 숫자는 아랫첨자로 변환
def format_name(formula):
    # + 또는 - 개수 세기
    plus_count = formula.count('+')
    minus_count = formula.count('-')
    
    # + 또는 -를 제거한 base formula
    base_formula = formula.replace('+', '').replace('-', '')
    
    # 숫자를 아랫첨자로 변환
    formatted_formula = ''
    for char in base_formula:
        if char.isdigit():
            formatted_formula += SUBSCRIPT_NUMS[char]
        else:
            formatted_formula += char
    
    # 전하 표시 추가
    if plus_count > 0:
        if plus_count == 1:
            return formatted_formula + '⁺'
        else:
            # 전하 숫자를 윗첨자로 변환
            superscript_num = SUPERSCRIPT_NUMS.get(str(plus_count), str(plus_count))
            return formatted_formula + superscript_num + '⁺'
    elif minus_count > 0:
        if minus_count == 1:
            return formatted_formula + '⁻'
        else:
            # 전하 숫자를 윗첨자로 변환
            superscript_num = SUPERSCRIPT_NUMS.get(str(minus_count), str(minus_count))
            return formatted_formula + superscript_num + '⁻'
    else:
        return formatted_formula

# conc 컬럼만 과학적 표기법으로 formatting하는 함수
def format_df_for_display(df):
    df_copy = df.copy()
    
    # 숫자 컬럼들을 소숫점 2자리로 formatting
    numeric_cols = df_copy.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        if col == 'conc':
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.0e}")
        else:
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.2f}")
    
    return df_copy

if __name__ == "__main__":
    # units
    kjmol = 96.485
    calmol = 23.061

    # constants
    kb = 8.617e-5 # eV/K
    T = 298.15 # K
    const = kb * T * np.log(10) # 0.0592 eV
    water = 56.690/calmol

    # gas
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

    gh = gh2 / 2
    go = gh2o - gh2
    goh = gh2o - gh2 / 2
    gooh = 2 * gh2o - 1.5 * gh2

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

    main()