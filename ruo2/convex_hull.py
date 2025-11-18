#!/usr/bin/env python3
"""
Convex Hull Analysis for Phase Diagrams

이 스크립트는 /Users/jiuy97/Desktop/3_RuO2/4_high_valence 폴더의 데이터를 읽어서
convex hull을 계산하고 분석합니다.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from ase.io import read
from ase import Atoms
from collections import defaultdict
import warnings
import re

warnings.filterwarnings('ignore')

# 금속 인덱스 매핑
METAL_MAP = {'Ru': 0, 'Hf': 1, 'Ta': 2, 'W': 3, 'Re': 4, 'Os': 5}


def load_compounds_from_folder(folder_path: Path, folder_name: str) -> List[Dict]:
    """폴더에서 화합물 데이터를 로드합니다."""
    compounds = []
    if not folder_path.exists():
        return compounds
    
    for subfolder in sorted(folder_path.iterdir()):
        if subfolder.is_dir():
            json_path = subfolder / "final_with_calculator.json"
            if json_path.exists():
                try:
                    atoms = read(str(json_path))
                    energy = atoms.get_total_energy()
                    composition = defaultdict(int)
                    for symbol in atoms.get_chemical_symbols():
                        composition[symbol] += 1
                    
                    compounds.append({
                        'name': subfolder.name,
                        'composition': dict(composition),
                        'energy': energy,
                        'folder': folder_name
                    })
                    print(f"Loaded: {subfolder.name} - Energy: {energy:.6f} eV")
                except Exception as e:
                    print(f"Error reading {json_path}: {e}")
    
    return compounds


def load_high_valence_data(base_path: str = "/Users/jiuy97/Desktop/3_RuO2/4_high_valence"):
    """high_valence 폴더에서 모든 데이터를 로드합니다."""
    base = Path(base_path)
    compounds = []
    
    for folder_name in ["1_M-RuO2", "2_MO2", "3_MxOy"]:
        folder_path = base / folder_name
        compounds.extend(load_compounds_from_folder(folder_path, folder_name))
    
    return compounds


def filter_compounds_for_mo2_pair(compounds: List[Dict], mo2_pair: Tuple[str, str]) -> List[Dict]:
    """특정 MO2 쌍에 관련된 화합물만 필터링합니다."""
    m1, m2 = mo2_pair
    m1_idx = METAL_MAP.get(m1)
    m2_idx = METAL_MAP.get(m2)
    
    if m1_idx is None or m2_idx is None:
        print(f"경고: 알 수 없는 금속: {m1} 또는 {m2}")
        return []
    
    filtered = []
    for compound in compounds:
        folder = compound['folder']
        try:
            folder_idx = int(compound['name'].split('_')[0])
        except (ValueError, IndexError):
            continue
        
        if folder == '2_MO2' and (folder_idx == m1_idx or folder_idx == m2_idx):
            filtered.append(compound)
        elif folder == '1_M-RuO2':
            if (m1 == 'Ru' and folder_idx == m2_idx) or (m2 == 'Ru' and folder_idx == m1_idx):
                filtered.append(compound)
    
    return filtered


def calculate_formation_energy_with_mo2_reference(compounds: List[Dict], 
                                                   m1: str, m2: str) -> List[Dict]:
    """MO2 쌍을 reference로 사용하여 형성에너지를 계산합니다."""
    # M1O2와 M2O2 찾기
    m1o2_energy = None
    m2o2_energy = None
    m1o2_atoms = None
    m2o2_atoms = None
    
    # 디버깅: 모든 화합물 출력
    print(f"\n화합물 목록:")
    for compound in compounds:
        comp = compound['composition']
        print(f"  {compound['name']}: {comp}, Energy: {compound['energy']:.6f} eV")
    
    for compound in compounds:
        comp = compound['composition']
        # 순수 MO2인지 확인 (해당 금속과 O만 포함, 다른 원소 없음)
        other_elements = [el for el in comp.keys() if el not in [m1, m2, 'O']]
        if len(other_elements) == 0 and 'O' in comp:
            # M1O2 찾기
            if m1 in comp and m2 not in comp:
                # O와 M1만 있고, O:M1 비율이 약 2:1인지 확인
                o_count = comp.get('O', 0)
                m1_count = comp.get(m1, 0)
                if m1_count > 0 and abs(o_count / m1_count - 2.0) < 0.5:  # 유연한 조건
                    m1o2_energy = compound['energy']
                    m1o2_atoms = sum(comp.values())
                    print(f"찾음: {m1}O2 = {compound['name']}, {comp}")
            # M2O2 찾기
            elif m2 in comp and m1 not in comp:
                o_count = comp.get('O', 0)
                m2_count = comp.get(m2, 0)
                if m2_count > 0 and abs(o_count / m2_count - 2.0) < 0.5:
                    m2o2_energy = compound['energy']
                    m2o2_atoms = sum(comp.values())
                    print(f"찾음: {m2}O2 = {compound['name']}, {comp}")
    
    if m1o2_energy is None or m2o2_energy is None:
        print(f"\n경고: {m1}O2 또는 {m2}O2를 찾을 수 없습니다.")
        print(f"  {m1}O2 찾음: {m1o2_energy is not None}")
        print(f"  {m2}O2 찾음: {m2o2_energy is not None}")
        return compounds
    
    print(f"\nReference: {m1}O2 = {m1o2_energy:.6f} eV ({m1o2_atoms} atoms)")
    print(f"Reference: {m2}O2 = {m2o2_energy:.6f} eV ({m2o2_atoms} atoms)")
    
    # 각 화합물의 형성에너지 계산
    for compound in compounds:
        comp = compound['composition']
        total_atoms = sum(comp.values())
        
        # M2의 조성 비율 계산
        m2_fraction = comp.get(m2, 0) / total_atoms if m2 in comp else 0.0
        
        # Reference 에너지 계산 (선형 보간)
        # 예: Ru_x Hf_(1-x) O_y 조성의 경우
        # reference = x * E(RuO2) + (1-x) * E(HfO2) (원자 수에 맞게 조정)
        # 하지만 더 정확하게는 조성에 따라 계산해야 함
        
        # 간단한 방법: M2의 원자 비율에 따라 reference 에너지 계산
        m1_fraction = comp.get(m1, 0) / total_atoms if m1 in comp else 0.0
        total_metal_fraction = m1_fraction + m2_fraction
        
        if total_metal_fraction > 0:
            # 금속 비율 정규화
            norm_m1_frac = m1_fraction / total_metal_fraction
            norm_m2_frac = m2_fraction / total_metal_fraction
            
            # Reference 에너지 (원자 수에 맞게 스케일링)
            ref_energy = (norm_m1_frac * m1o2_energy / m1o2_atoms + 
                         norm_m2_frac * m2o2_energy / m2o2_atoms) * total_atoms
        else:
            ref_energy = 0.0
        
        # 형성에너지 = (총에너지 - reference에너지) / 원자 수
        formation_energy = (compound['energy'] - ref_energy) / total_atoms
        compound['formation_energy'] = formation_energy
    
    return compounds


def plot_binary_phase_diagram(compounds: List[Dict], m1: str, m2: str,
                              output_file: Optional[str], show_plot: bool):
    """2원계 phase diagram을 그립니다."""
    # 조성과 형성에너지 추출
    compositions, formation_energies, labels = [], [], []
    
    print(f"\n플로팅 데이터:")
    for compound in compounds:
        comp = compound['composition']
        total_atoms = sum(comp.values())
        
        # M2의 조성 비율 계산
        fraction = comp.get(m2, 0) / total_atoms if m2 in comp else 0.0
        formation_energy = compound.get('formation_energy', 0.0)
        
        compositions.append(fraction)
        formation_energies.append(formation_energy)
        
        # 라벨 생성
        comp_str = ''.join([f"{el}{count}" if count > 1 else el 
                           for el, count in sorted(comp.items())])
        labels.append(comp_str)
        print(f"  {comp_str}: x={fraction:.3f}, E_form={formation_energy:.6f} eV/atom")
    
    if len(compositions) == 0:
        print("경고: 플로팅할 데이터가 없습니다.")
        return
    
    # Convex hull 계산
    points = np.column_stack([compositions, formation_energies])
    print(f"\nPoints shape: {points.shape}")
    print(f"Points:\n{points}")
    
    try:
        hull = ConvexHull(points)
        print(f"Convex hull 계산 성공: {len(hull.vertices)} vertices")
    except Exception as e:
        print(f"Convex hull 계산 실패: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # 플롯
    fig, ax = plt.subplots(figsize=(10, 8))
    stable_indices = set(hull.vertices)
    
    # 안정적인 상과 불안정한 상 구분하여 플롯
    stable_x = [compositions[i] for i in stable_indices]
    stable_y = [formation_energies[i] for i in stable_indices]
    unstable_x = [compositions[i] for i in range(len(compositions)) if i not in stable_indices]
    unstable_y = [formation_energies[i] for i in range(len(compositions)) if i not in stable_indices]
    
    if stable_x:
        ax.scatter(stable_x, stable_y, marker='o', s=150, edgecolor='black', 
                  facecolor='green', zorder=3, label='Stable')
    if unstable_x:
        ax.scatter(unstable_x, unstable_y, marker='D', s=100, edgecolor='black', 
                  facecolor='red', alpha=0.6, zorder=3, label='Above hull')
    
    # Convex hull 그리기
    hull_vertices = points[hull.vertices]
    hull_vertices_sorted = hull_vertices[hull_vertices[:, 0].argsort()]
    hull_vertices_closed = np.vstack([hull_vertices_sorted, hull_vertices_sorted[0]])
    ax.plot(hull_vertices_closed[:, 0], hull_vertices_closed[:, 1],
           'k-', linewidth=2, zorder=2, label='Convex hull')
    
    # 라벨 추가
    for x, y, label in zip(compositions, formation_energies, labels):
        ax.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points', 
                   fontsize=9, ha='left')
    
    ax.set_xlabel(f'Composition Fraction ({m2})', fontsize=12)
    ax.set_ylabel('Formation Energy (eV/atom)', fontsize=12)
    ax.set_title(f'{m1}O2 - {m2}O2 Phase Diagram', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Phase diagram이 {output_file}에 저장되었습니다.")
    
    if show_plot:
        plt.show()
    
    # 안정적인 상 출력
    print(f"\n=== 안정적인 상 (Stable Phases) ===")
    for i, idx in enumerate(sorted(stable_indices), 1):
        print(f"{i}. {labels[idx]}: {formation_energies[idx]:.6f} eV/atom")


def plot_combined_mo2_diagram(base_path: str, pairs: List[Tuple[str, str]],
                              output_file: str, show_plot: bool = False):
    """5개 MO2 쌍을 하나의 그래프에 합쳐서 그립니다."""
    print("=== Combined MO2 Pairs Convex Hull Analysis ===")
    print(f"데이터 경로: {base_path}\n")
    
    # 데이터 로드
    all_compounds = load_high_valence_data(base_path)
    print(f"총 {len(all_compounds)}개의 화합물을 로드했습니다.\n")
    
    if len(all_compounds) == 0:
        print("로드된 화합물이 없습니다.")
        return
    
    # 색상과 마커 설정
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # 파, 주황, 초록, 빨강, 보라
    markers = ['o', 's', '^', 'D', 'v']
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 각 쌍에 대해 처리
    for idx, (m1, m2) in enumerate(pairs):
        print(f"{'='*60}")
        print(f"Processing: {m1}O2 - {m2}O2")
        print(f"{'='*60}\n")
        
        compounds = filter_compounds_for_mo2_pair(all_compounds, (m1, m2))
        print(f"{m1}O2-{m2}O2 시스템: {len(compounds)}개의 화합물 발견")
        
        if len(compounds) < 2:
            print(f"경고: {m1}O2-{m2}O2 시스템에 충분한 데이터가 없습니다.\n")
            continue
        
        # 형성에너지 계산 (MO2 쌍을 reference로 사용)
        compounds = calculate_formation_energy_with_mo2_reference(compounds, m1, m2)
        
        # 조성과 형성에너지 추출
        compositions, formation_energies, labels = [], [], []
        
        for compound in compounds:
            comp = compound['composition']
            total_atoms = sum(comp.values())
            fraction = comp.get(m2, 0) / total_atoms if m2 in comp else 0.0
            compositions.append(fraction)
            formation_energies.append(compound.get('formation_energy', 0.0))
            
            comp_str = ''.join([f"{el}{count}" if count > 1 else el 
                               for el, count in sorted(comp.items())])
            labels.append(comp_str)
        
        if len(compositions) == 0:
            continue
        
        # Convex hull 계산
        points = np.column_stack([compositions, formation_energies])
        try:
            hull = ConvexHull(points)
            stable_indices = set(hull.vertices)
            
            # 안정적인 상과 불안정한 상 구분
            stable_x = [compositions[i] for i in stable_indices]
            stable_y = [formation_energies[i] for i in stable_indices]
            unstable_x = [compositions[i] for i in range(len(compositions)) if i not in stable_indices]
            unstable_y = [formation_energies[i] for i in range(len(compositions)) if i not in stable_indices]
            
            # 안정적인 상 플롯
            if stable_x:
                ax.scatter(stable_x, stable_y, marker=markers[idx], s=100, 
                          edgecolor='black', facecolor=colors[idx], 
                          zorder=3, label=f'{m2}', alpha=0.8)
            
            # 불안정한 상 플롯
            if unstable_x:
                ax.scatter(unstable_x, unstable_y, marker=markers[idx], s=60, 
                          edgecolor='black', facecolor=colors[idx], 
                          alpha=0.4, zorder=2)
            
            # Convex hull 그리기
            hull_vertices = points[hull.vertices]
            hull_vertices_sorted = hull_vertices[hull_vertices[:, 0].argsort()]
            hull_vertices_closed = np.vstack([hull_vertices_sorted, hull_vertices_sorted[0]])
            ax.plot(hull_vertices_closed[:, 0], hull_vertices_closed[:, 1],
                   color=colors[idx], linewidth=2, zorder=1, alpha=0.7, linestyle='--')
            
        except Exception as e:
            print(f"Convex hull 계산 실패: {e}")
            continue
        
        print()
    
    # 축 설정
    ax.set_xlabel('Composition Fraction (M)', fontsize=14)
    ax.set_ylabel('Formation Energy (eV/atom)', fontsize=14)
    ax.set_title('RuO2 - MO2 Phase Diagrams (Combined)', fontsize=16)
    ax.legend(loc='best', fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Combined phase diagram이 {output_file}에 저장되었습니다.")
    
    if show_plot:
        plt.show()


def plot_mo2_pairs(base_path: str, pairs: List[Tuple[str, str]], 
                   output_dir: str = ".", show_plot: bool = False):
    """여러 MO2 쌍에 대한 phase diagram을 생성합니다."""
    print("=== MO2 Pairs Convex Hull Analysis ===")
    print(f"데이터 경로: {base_path}\n")
    
    # 데이터 로드
    all_compounds = load_high_valence_data(base_path)
    print(f"총 {len(all_compounds)}개의 화합물을 로드했습니다.\n")
    
    if len(all_compounds) == 0:
        print("로드된 화합물이 없습니다.")
        return
    
    # 각 쌍에 대해 처리
    for m1, m2 in pairs:
        print(f"{'='*60}")
        print(f"Processing: {m1}O2 - {m2}O2")
        print(f"{'='*60}\n")
        
        compounds = filter_compounds_for_mo2_pair(all_compounds, (m1, m2))
        print(f"{m1}O2-{m2}O2 시스템: {len(compounds)}개의 화합물 발견")
        
        if len(compounds) < 2:
            print(f"경고: {m1}O2-{m2}O2 시스템에 충분한 데이터가 없습니다.\n")
            continue
        
        # 형성에너지 계산 (MO2 쌍을 reference로 사용)
        compounds = calculate_formation_energy_with_mo2_reference(compounds, m1, m2)
        
        # Phase diagram 그리기
        output_file = f"{output_dir}/{m1}O2-{m2}O2_convex_hull.png"
        plot_binary_phase_diagram(compounds, m1, m2, output_file, show_plot)
        print()
    
    # Combined 그래프 생성
    print(f"\n{'='*60}")
    print("Creating combined diagram...")
    print(f"{'='*60}\n")
    combined_output = f"{output_dir}/RuO2-MO2_combined_convex_hull.png"
    plot_combined_mo2_diagram(base_path, pairs, combined_output, show_plot)


def main():
    """메인 함수"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='High Valence Compounds Convex Hull Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예제:
  # MO2 쌍들 분석 (RuO2-HfO2, RuO2-TaO2, 등)
  python convex_hull.py --mo2-pairs
  
  # 출력 디렉토리 지정
  python convex_hull.py --mo2-pairs --output-dir ./results
        """
    )
    
    parser.add_argument('--base-path', type=str,
                       default="/Users/jiuy97/Desktop/3_RuO2/4_high_valence",
                       help='데이터 폴더 경로')
    parser.add_argument('--output-dir', type=str, default='.',
                       help='출력 디렉토리')
    parser.add_argument('--no-show', action='store_true',
                       help='플롯을 화면에 표시하지 않음')
    parser.add_argument('--mo2-pairs', action='store_true',
                       help='MO2 쌍들에 대한 phase diagram 생성')
    
    args = parser.parse_args()
    
    # MO2 쌍 분석 모드
    pairs = [
        ('Ru', 'Hf'), ('Ru', 'Ta'), ('Ru', 'W'), 
        ('Ru', 'Re'), ('Ru', 'Os'),
    ]
    plot_mo2_pairs(args.base_path, pairs, args.output_dir, not args.no_show)
    print("\n모든 MO2 쌍 분석 완료!")
    print("개별 그래프와 combined 그래프가 생성되었습니다.")


if __name__ == "__main__":
    main()
