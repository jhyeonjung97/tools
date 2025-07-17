#!/bin/env python

import os
import socket
import numpy as np
import pandas as pd
from ase.io import read, write
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 서버 주소 가져오기
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/8_V_slab'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/8_V_slab'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/8_V_slab'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")
save_path = os.path.join(root, 'figures')
os.makedirs(save_path, exist_ok=True)

# 참조 에너지 값들
h2o = -14.23919983
h2 = -6.77409008

zpeh2o = 0.558
zpeh2 = 0.273

cph2o = 0.10
cph2 = 0.09

gh2 = h2 + zpeh2 + cph2
gh2o = h2o + zpeh2o + cph2o
go = gh2o - gh2  # O 원자의 Gibbs 자유에너지
goh = gh2o - gh2/2  # OH의 Gibbs 자유에너지

coords_data = [
    {'coord': 'WZ', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.667, 'coord_dir': '1_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
    {'coord': 'ZB', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.833, 'coord_dir': '2_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
    {'coord': 'TN', 'CN': 4, 'OS': 2, 'MN': 4, 'surfox': 1.933, 'coord_dir': '3_SquarePlanar_TN', 'zorder': 3, 'marker': 'o', 'color': 'dodgerblue',},
    {'coord': 'PD', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 1.667, 'coord_dir': '4_SquarePlanar_PD', 'zorder': 2, 'marker': 'o', 'color': 'deepskyblue',},
    {'coord': 'NB', 'CN': 4, 'OS': 2, 'MN': 6, 'surfox': 2.667, 'coord_dir': '5_SquarePlanar_NB', 'zorder': 1, 'marker': 's', 'color': 'limegreen',},
    {'coord': 'RS', 'CN': 6, 'OS': 2, 'MN': 2, 'surfox': 1.933, 'coord_dir': '6_Octahedral_RS',   'zorder': 6, 'marker': 'd', 'color': 'orchid',},
    {'coord': 'LT', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '7_Pyramidal_LT',    'zorder': 0, 'marker': 'h', 'color': 'silver',},
    {'coord': 'WW', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '8_Tetrahedral_WZ',  'zorder': 5, 'marker': '>', 'color': 'darkorange',},
    {'coord': 'ZZ', 'CN': 4, 'OS': 2, 'MN': 2, 'surfox': 2.000, 'coord_dir': '9_Tetrahedral_ZB',  'zorder': 4, 'marker': '<', 'color': 'gold',},
]

coords = pd.DataFrame(coords_data).set_index('coord')
coords.index.name = None

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
    'fm': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
}

def get_energy(atoms_path):
    """에너지를 계산합니다."""
    atoms = read(atoms_path)
    return atoms.get_total_energy()

def get_clean_surface_energy(path):
    """깨끗한 표면의 에너지를 계산하고 x축 뷰 이미지를 저장합니다."""
    clean_path = os.path.join(path, 'clean', 'final_with_calculator.json')
    if os.path.exists(clean_path):
        atoms = read(clean_path)
        atoms_for_image = atoms.copy()
        atoms_for_image = atoms_for_image.repeat((2, 2, 1))
        slab_path = os.path.join(path, 'slab.png')
        write(slab_path, atoms_for_image, rotation='-90x', show_unit_cell=2)
        return atoms.get_total_energy()
    return None

def get_calculation_status(path):
    """계산 폴더의 상태를 확인합니다."""
    if not os.path.exists(os.path.join(path, 'submit.sh')):
        return '-'
    if os.path.exists(os.path.join(path, 'unmatched')):
        return 'unmtch'
    elif os.path.exists(os.path.join(path, 'unstable')):
        return 'unstbl'
    elif os.path.exists(os.path.join(path, 'DONE')):
        return 'done'
    else:
        return 'submit'

def get_lowest_energy_adsorption(adsorption_paths, clean_energy):
    """여러 흡착 구조 중 가장 낮은 에너지를 가진 것을 선택하고 흡착 에너지를 계산합니다."""
    energies = {}
    for path in adsorption_paths:
        # 흡착물 폴더에 DONE 파일이 있는지 확인하고 unmatched 파일이 없는지 확인
        ads_dir = os.path.dirname(path)
        if not os.path.exists(os.path.join(ads_dir, 'DONE')) or os.path.exists(os.path.join(ads_dir, 'unmatched')):
            continue
            
        energy = get_energy(path)
        # o1, o2와 같은 숫자를 제거하여 기본 흡착물 종류를 얻습니다
        base_ads = ''.join([c for c in os.path.basename(ads_dir) if not c.isdigit()])
        
        # 흡착 에너지 계산
        if base_ads == 'o':
            ads_energy = energy - clean_energy - go + 0.28931784
        elif base_ads == 'oh':
            ads_energy = energy - clean_energy - goh + 0.46951875
        else:
            continue
            
        if base_ads not in energies or ads_energy < energies[base_ads]['energy']:
            energies[base_ads] = {'energy': ads_energy, 'path': path}
    return energies

def create_progress_file():
    """계산 진행 상황을 추적하는 TSV 파일을 생성합니다."""
    # 컬럼 생성
    columns = ['coord', 'row', 'ads']
    columns.extend([str(i).zfill(2) for i in range(13)])  # 00부터 12까지
    
    # 데이터프레임 생성
    df = pd.DataFrame(columns=columns)
    
    # 각 coord와 row에 대해
    for coord in coords.index:
        coord_dir = coords.loc[coord, 'coord_dir']
        for row in ['3d', '4d', '5d']:
            # print(f"\nProcessing {coord} {row}...")
            # clean 상태 추가
            clean_status = []
            for i in range(13):
                path = os.path.join(root, coord_dir, row, f'{i:02d}_{metals[row][i]}', 'clean')
                status = get_calculation_status(path)
                clean_status.append(status)
            df.loc[len(df)] = [coord, row, 'clean'] + clean_status
            
            # o, o1, o2 상태 추가
            for ads in ['o', 'o1', 'o2']:
                ads_status = []
                for i in range(13):
                    path = os.path.join(root, coord_dir, row, f'{i:02d}_{metals[row][i]}', ads)
                    status = get_calculation_status(path)
                    ads_status.append(status)
                df.loc[len(df)] = [coord, row, ads] + ads_status
            
            # oh, oh1, oh2 상태 추가
            for ads in ['oh', 'oh1', 'oh2']:
                ads_status = []
                for i in range(13):
                    path = os.path.join(root, coord_dir, row, f'{i:02d}_{metals[row][i]}', ads)
                    status = get_calculation_status(path)
                    ads_status.append(status)
                df.loc[len(df)] = [coord, row, ads] + ads_status
    
    # 모든 금속 데이터가 비어있는 행 삭제
    metal_columns = [str(i).zfill(2) for i in range(13)]
    df = df[df[metal_columns].apply(lambda x: x != '-').any(axis=1)]
    
    # TSV 파일로 저장
    output_path = os.path.join(save_path, 'calculation_progress.tsv')
    df.to_csv(output_path, sep='\t', index=False)
    # print(f"\nProgress file saved to: {output_path}")
    return df

def main():
    df = pd.DataFrame(columns=['coord', 'row', 'numb', 'metal', 'surfox', 'o_energy', 'oh_energy'])
    
    # 모든 가능한 경로 패턴을 찾습니다
    pattern = os.path.join(root, '*_*_*', '*', '*_*')
    all_paths = glob.glob(pattern)
    
    for path in all_paths:
        # DONE 파일이 있는지 확인
        if not os.path.exists(os.path.join(path, 'clean', 'DONE')) or os.path.exists(os.path.join(path, 'clean', 'unmatched')):
            continue
            
        # 경로에서 정보 추출
        parts = path.split('/')
        coord_dir = parts[-3]  # *_*_* 형식
        row = parts[-2]        # row (3d, 4d, 5d)
        metal_dir = parts[-1]  # *_* 형식
        
        # coord 정보 추출
        coord = None
        for c in coords.index:
            if c in coord_dir:
                coord = c
                break
        if coord is None:
            continue
            
        # metal 정보 추출
        metal = None
        for m in metals[row]:
            if m in metal_dir:
                metal = m
                break
        if metal is None:
            continue
            
        # metal의 index 추출
        numb = str(metals[row].index(metal)).zfill(2)
        item = coord + row + numb
        surfox = coords.loc[coord, 'surfox']
        
        # 깨끗한 표면 에너지 계산
        clean_energy = get_clean_surface_energy(path)
        if clean_energy is None:
            continue
        
        # 흡착 구조 찾기
        adsorption_paths = glob.glob(os.path.join(path, '*', 'final_with_calculator.json'))
        if not adsorption_paths:
            continue
            
        # 가장 낮은 에너지를 가진 흡착 구조 선택
        energies = get_lowest_energy_adsorption(adsorption_paths, clean_energy)
        
        # 데이터프레임에 추가
        df.loc[item, ['coord', 'row', 'numb', 'metal', 'surfox']] = coord, row, numb, metal, surfox
        if 'o' in energies:
            df.loc[item, 'o_energy'] = energies['o']['energy'] + 0.2975
        if 'oh' in energies:
            df.loc[item, 'oh_energy'] = energies['oh']['energy'] - 0.014
    
    # 결과 저장
    df = df.dropna(subset=['o_energy', 'oh_energy'], how='all')
    df.to_csv(f'{save_path}/slab_data.csv', sep=',')
    df[['o_energy', 'oh_energy']] = df[['o_energy', 'oh_energy']].astype(float).round(2)
    df.to_csv(f'{save_path}/slab_data.tsv', sep='\t', float_format='%.2f')
    
    # coord, row, numb 순으로 정렬 (numb를 int로 변환)
    df['numb'] = df['numb'].astype(int)
    df = df.sort_values(['coord', 'row', 'numb'])
    print(df)
    
    # 계산 진행 상황 파일 생성
    # print("\nCreating calculation progress file...")
    progress_df = create_progress_file()
    # print("\nCalculation Progress:")
    print(progress_df)
    
    # 그래프 그리기
    # plot_by_metal_row(df, save_path)
    # plot_by_coordination(df, save_path)

def plot_by_metal_row(df, save_path):
    """각 금속 행(row)별로 흡착 에너지를 그래프로 그립니다."""
    for row in ['fm', '3d', '4d', '5d']:
        metals_list = metals[row]
        for ads in ['o_energy', 'oh_energy']:
            plt.figure(figsize=(4, 3))
            for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT', 'WW', 'ZZ']:
                zorder = coords.loc[coord, 'zorder']
                marker = coords.loc[coord, 'marker']
                color = coords.loc[coord, 'color']
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                # metals_list 순서에 맞게 reindex
                subset = subset.set_index('metal').reindex(metals_list).reset_index()
                plt.plot(np.arange(len(metals_list)), subset[ads], marker=marker, color=color, 
                        linestyle='-', label=coord, zorder=zorder)
            plt.xticks(np.arange(len(metals_list)), metals_list)
            plt.xlabel("Metal Index")
            plt.ylabel(f"{ads.split('_')[0].upper()} Adsorption Energy (eV)")
            plt.xlim(-0.5, 12.5)
            plt.legend()
            plt.tight_layout()
            png_name = f"slab_{row}_{ads}.png"
            plt.savefig(f"{save_path}/{png_name}", transparent=True, dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")

def plot_by_coordination(df, save_path):
    """각 좌표계(coord)별로 흡착 에너지를 그래프로 그립니다."""
    for coord in ['WZ', 'ZB', 'TN', 'PD', 'NB', 'RS', 'LT', 'WW', 'ZZ']:
        for ads in ['o_energy', 'oh_energy']:
            marker = coords.loc[coord, 'marker']
            base_color = coords.loc[coord, 'color']
            cmap = mcolors.LinearSegmentedColormap.from_list(f'cmap_{base_color}', [base_color, 'white'])
            colors = cmap(np.linspace(0.0, 0.6, 3))
            plt.figure(figsize=(4, 3))
            for r, row in enumerate(['3d', '4d', '5d']):
                metals_list = metals[row]
                color = 'lightgray' if row == 'fm' else colors[r]
                subset = df[(df['coord'] == coord) & (df['row'] == row)]
                subset = subset.set_index('metal').reindex(metals_list).reset_index()
                plt.plot(np.arange(len(metals_list)), subset[ads], marker=marker, color=color, 
                        linestyle='-', label=row)
            plt.xticks(np.arange(len(metals['3d'])), [f'{a}\n{b}\n{c}' for a, b, c in zip(metals['3d'], metals['4d'], metals['5d'])])
            plt.xlabel("Metal Index")
            plt.ylabel(f"{ads.split('_')[0].upper()} Adsorption Energy (eV)")
            plt.xlim(-0.5, 12.5)
            plt.legend()
            plt.tight_layout()
            png_name = f"slab_{coord}_{ads}.png"
            plt.savefig(f"{save_path}/{png_name}", transparent=True, dpi=300)
            plt.close()
            print(f"Figure saved as {png_name}")

if __name__ == '__main__':
    main() 