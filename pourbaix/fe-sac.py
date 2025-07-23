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

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
parser.add_argument('--json-dir', type=str, default='.', help='json 파일이 있는 폴더 경로 (기본값: 현재 폴더)')
parser.add_argument('--csv-dir', type=str, default='.', help='label.csv가 있는 폴더 경로 (기본값: 현재 폴더)')
parser.add_argument('--conc', type=float, default=10**-6, help='concentration (기본값: 10^-6)')
parser.add_argument('--oh', action='store_true', help='label.csv의 세 번째 컬럼을 OH 개수로 사용')
parser.add_argument('--bulk', action='store_true', help='bulk 모드')
parser.add_argument('--ads', type=str, nargs='*', default=[], help='adsorbate 원소들 (예: --ads N O)')
parser.add_argument('--ph', type=int, default=0, help='pH value for the plot (default: 0)')
parser.add_argument('--tick', type=float, default=0.01, help='Tick size for pH and U ranges (default: 0.01)')
parser.add_argument('--figx', type=float, default=4, help='Figure width in inches (default: 6)')
parser.add_argument('--figy', type=float, default=3, help='Figure height in inches (default: 7)')
parser.add_argument('--show', action='store_true', help='Show the plot')
args = parser.parse_args()

is_bulk = args.bulk
ads_elements = args.ads

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
        row['name'] = file_labels.get(os.path.basename(json_file), json_file)
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
        oh_count = 0
        if args.oh and os.path.exists(label_csv_path):
            surf_name = surfs[k]['name']
            label_df = pd.read_csv(label_csv_path, header=None)
            matching_row = label_df[label_df[1] == surf_name]  # 두 번째 컬럼(인덱스 1)과 매칭
            if not matching_row.empty and len(matching_row.iloc[0]) > 2:
                oh_count = float(matching_row.iloc[0][2])  # 세 번째 컬럼(인덱스 2)에서 가져오기
                # print(f"Found OH count: {oh_count}")
            else:
                print(f"No match found for '{surf_name}' or insufficient columns")        
        formation_energy_correction = (
            - (surfs[k]['H'] - oh_count) * (gh - dgh)
            - (surfs[k]['O'] - oh_count) * (go - dgo)
            - oh_count * (goh - dgoh)
        )
        surfs[k]['E_DFT'] = surfs[k]['E_DFT'] - reference_surface_energy + formation_energy_correction
    
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
                        row['name'] = ion_formula
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            row['name'] = f'1/{int(el_count)} ' + row['name']
                        
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
                        row['name'] = solid_formula
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            row['name'] = f'1/{int(el_count)} ' + row['name']
                        
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
                        row['name'] = gas_formula
                        
                        # el 개수로 normalization
                        el_count = reduced_dict.get(el, 1)  # el이 없으면 1로 설정
                        if el_count > 1:
                            row['E_DFT'] = row['E_DFT'] / el_count
                            row['e'] = row['e'] / el_count
                            for elem in remaining_elements:
                                row[elem] = row[elem] / el_count
                            row['name'] = f'1/{int(el_count)} ' + row['name']
                        
                        gases.append(row)
                    except:
                        print(f"    {gas_formula}: parsing failed, energy: {energy}")
        else:
            print(f"\n{el}: No thermodynamic data found")

    nsurfs, nions, nsolids, ngases = len(surfs), len(ions), len(solids), len(gases)
    surfs_df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    ions_df = pd.DataFrame(ions, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    solids_df = pd.DataFrame(solids, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    gases_df = pd.DataFrame(gases, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    
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
                                else:
                                    new_surf[key] = surfs[k][key] + gas[key]
                            new_surfs.append(new_surf)

    # new_surfs를 surfs에 추가
    surfs.extend(new_surfs)
    # unique_elements가 모두 0이 아닌 것들만 남기기
    surfs = [surf for surf in surfs if not all(surf[elem] == 0 for elem in unique_elements)]
    nsurfs = len(surfs)
    df = pd.DataFrame(surfs, columns=['E_DFT', 'e'] + remaining_elements + ['conc', 'name'])
    print(format_df_for_display(df))
    print(f"After adding combinations: {nsurfs} surfs entries")

    # pH and potential grid 설정
    tick = args.tick if hasattr(args, 'tick') else 0.01
    pHmin, pHmax = -2, 16
    pHrange = np.arange(pHmin, pHmax + tick, tick)
    Umin, Umax = -1.0, 3.0
    Urange = np.arange(Umin, Umax + 0.06 * 14, tick)
    target_pH = args.ph if hasattr(args, 'ph') else 0

    # dg 함수 정의 (Gibbs free energy 계산)
    def dg(k, pH, U, n_ref):
        surface_term = surfs[k]['E_DFT'] - surfs[n_ref]['E_DFT']
        U_coeff = surfs[k]['H'] - 2*surfs[k]['O'] - surfs[k]['e']
        pH_coeff = surfs[k]['H'] - 2*surfs[k]['O'] - surfs[k]['e']
        dg_value = surface_term + U_coeff*U + const*pH_coeff*pH + const*log10(surfs[k]['conc'])
        return dg_value

    # reference surface index 찾기
    n_ref = ref_surf_idx

    # lowest surfaces 계산
    lowest_surfaces = np.full((len(Urange), len(pHrange)), np.nan)

    pHindex = 0
    for pH in pHrange:
        Uindex = 0
        for U in Urange:
            values = []
            for k in range(nsurfs):
                value = dg(k, pH, U, n_ref=n_ref)
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
            value = dg(k, pH=target_pH, U=U, n_ref=n_ref)
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