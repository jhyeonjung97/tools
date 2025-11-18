#!/usr/bin/env python3
"""
Materials Project Phase Diagram Generator

이 스크립트는 Materials Project API를 사용하여 phase diagram을 생성합니다.
Materials Project 웹사이트에서 보는 것과 동일한 형식의 phase diagram을 그립니다.
"""

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from contextlib import nullcontext
import warnings

# MPRester import - 최신 버전과 구버전 모두 지원
try:
    from mp_api.client import MPRester
except ImportError:
    try:
        from pymatgen.ext.matproj import MPRester
    except ImportError:
        raise ImportError("MPRester를 찾을 수 없습니다. mp-api 또는 pymatgen을 설치하세요.")

warnings.filterwarnings('ignore')


def get_phase_diagram(elements, api_key=None, ref_material_ids=None):
    """
    Materials Project에서 phase diagram 데이터를 가져옵니다.
    
    Args:
        elements: 원소 리스트 (예: ['Re', 'O'] 또는 ['Re'])
        api_key: Materials Project API 키 (없으면 환경변수에서 가져옴)
        ref_material_ids: Reference state로 사용할 material ID 딕셔너리 (예: {'Re': 'mp-8', 'O': 'mp-12957'})
        
    Returns:
        PhaseDiagram 객체
    """
    if api_key is None:
        api_key = os.getenv('MAPI_KEY')
        if not api_key:
            sys.exit("Error: MAPI_KEY environment variable not set. "
                    "Please set your Materials Project API key as an environment variable.")
    
    print(f"Materials Project에서 {elements} 시스템의 데이터를 가져오는 중...")
    
    # 먼저 최신 mp_api를 시도하고, 실패하면 구버전 pymatgen MPRester 사용
    entries = None
    use_old_api = False
    
    try:
        with MPRester(api_key) as mpr:
            entries = mpr.get_entries_in_chemsys(elements)
            print(f"총 {len(entries)}개의 entries를 가져왔습니다.")
    except AttributeError as e:
        if "model_fields" in str(e):
            print("에러: mp_api 라이브러리 버전 호환성 문제가 발생했습니다.")
            print("pymatgen의 구버전 MPRester를 시도합니다...")
            use_old_api = True
        else:
            raise
    except Exception as e:
        print(f"에러 발생: {type(e).__name__}: {e}")
        print("pymatgen의 구버전 MPRester를 시도합니다...")
        use_old_api = True
    
    # 구버전 API 사용
    if use_old_api or entries is None:
        try:
            from pymatgen.ext.matproj import MPRester as OldMPRester
            with OldMPRester(api_key) as old_mpr:
                entries = old_mpr.get_entries_in_chemsys(elements)
                print(f"총 {len(entries)}개의 entries를 가져왔습니다 (구버전 API 사용).")
        except Exception as e2:
            print("구버전 API도 실패했습니다.")
            print("해결 방법:")
            print("1. mp-api를 업데이트: pip install --upgrade mp-api")
            print("2. 또는 다운그레이드: pip install 'mp-api<0.30.0'")
            print("3. 또는 pymatgen 설치: pip install pymatgen")
            sys.exit(f"구버전 API 에러: {e2}")
    
    # entries가 여전히 None이면 에러
    if entries is None:
        sys.exit("에러: entries를 가져올 수 없습니다.")
    
    # 나머지 코드는 entries만 사용하므로 mpr 객체가 필요 없음
    # nullcontext를 사용하여 with 블록 구조 유지
    with nullcontext():
        
        # 각 원소의 pure element entry 중 가장 낮은 에너지를 가진 것을 reference state로 선택
        best_ref_entries = {}
        for el in elements:
            # 해당 원소만 포함하는 모든 entry 찾기
            pure_entries = [e for e in entries 
                          if len(e.composition.elements) == 1 and 
                          e.composition.elements[0].symbol == el]
            
            if pure_entries:
                # energy_per_atom이 가장 낮은 entry 선택
                best_entry = min(pure_entries, key=lambda e: e.energy_per_atom)
                best_ref_entries[el] = best_entry
        
        # Reference state를 특정 material ID로 지정하는 경우
        if ref_material_ids:
            # 각 원소에 대해 지정된 material ID의 entry 찾기
            ref_entries_to_use = {}
            for el, mp_id in ref_material_ids.items():
                if el not in elements:
                    continue
                # material ID에서 숫자만 추출 (예: "mp-8" -> "8")
                mp_num = mp_id.replace('mp-', '').split('-')[0]
                
                # 해당 material ID의 entry 찾기
                found_entry = None
                for entry in entries:
                    entry_id = getattr(entry, 'entry_id', '')
                    material_id = getattr(entry, 'material_id', '')
                    
                    # entry_id에서 material ID 추출
                    if entry_id and 'mp-' in str(entry_id):
                        entry_mp_match = str(entry_id).split('-')
                        if len(entry_mp_match) >= 2:
                            entry_mp_num = entry_mp_match[1]
                            if entry_mp_num == mp_num:
                                found_entry = entry
                                break
                    
                    # material_id 속성에서 확인
                    if material_id and mp_num in str(material_id):
                        found_entry = entry
                        break
                
                if found_entry:
                    ref_entries_to_use[el] = found_entry
                    # 지정된 material ID를 best_ref_entries에 반영
                    best_ref_entries[el] = found_entry
        
        # Phase diagram 생성 전에 reference state를 강제로 설정
        # PhaseDiagram이 가장 낮은 에너지 entry를 선택하지 않을 수 있으므로,
        # entries를 수정하여 각 원소의 가장 낮은 에너지 entry만 남기고 나머지 pure element entries 제거
        entries_to_keep = []
        
        for entry in entries:
            # Pure element entry인지 확인
            is_pure_element = (len(entry.composition.elements) == 1 and 
                             entry.composition.elements[0].symbol in elements)
            
            if is_pure_element:
                el = entry.composition.elements[0].symbol
                # 가장 낮은 에너지 entry만 유지
                if el in best_ref_entries and entry == best_ref_entries[el]:
                    entries_to_keep.append(entry)
            else:
                # Pure element가 아닌 entry는 모두 유지
                entries_to_keep.append(entry)
        
        # 수정된 entries로 Phase diagram 생성
        pd = PhaseDiagram(entries_to_keep)
        print(f"안정적인 상: {len(pd.stable_entries)}개")
        
        return pd, entries_to_keep


def plot_phase_diagram(pd, elements, entries=None, output_file=None, show_plot=True, 
                      xlim=None, ylim=None, figsize=None, fontsize=14):
    """
    Phase diagram을 그립니다.
    
    Args:
        pd: PhaseDiagram 객체
        elements: 원소 리스트
        entries: 모든 entries 리스트 (불안정한 상 식별용, None이면 pd에서 추출 시도)
        output_file: 저장할 파일 경로 (None이면 저장 안함)
        show_plot: 플롯을 화면에 표시할지 여부
        xlim: x축 범위 튜플 (min, max), None이면 자동 설정
        ylim: y축 범위 튜플 (min, max), None이면 자동 설정
        figsize: 그림 사이즈 튜플 (width, height) 인치 단위, None이면 자동 설정
    """
    # 원소 개수에 따라 다른 plotter 사용
    num_elements = len(elements)
    
    if num_elements == 1:
        # 단일 원소: 1D plot (일반적으로 사용하지 않음)
        print("단일 원소 시스템입니다. 2원계 또는 3원계를 권장합니다.")
        plotter = PDPlotter(pd, show_unstable=True)
        plotter.get_plot()
        
    elif num_elements == 2:
        # 2원계: 2D plot - raw data로 직접 그리기
        print("2원계 phase diagram을 그리는 중...")
        
        # Figure 생성
        if figsize is not None:
            fig, ax = plt.subplots(figsize=figsize)
            print(f"그림 사이즈: {figsize[0]:.2f} x {figsize[1]:.2f} 인치")
        else:
            fig, ax = plt.subplots(figsize=(8, 6))
            figsize_current = fig.get_size_inches()
            print(f"그림 사이즈: {figsize_current[0]:.2f} x {figsize_current[1]:.2f} 인치")
        
        # 데이터 준비: 모든 entries에서 composition fraction과 formation energy 추출
        stable_x = []
        stable_y = []
        stable_labels = []
        unstable_x = []
        unstable_y = []
        unstable_labels = []
        unstable_e_above_hull = []  # energy_above_hull 값 저장
        
        # 모든 entries 확인
        all_entries_list = entries if entries is not None else pd.all_entries if hasattr(pd, 'all_entries') else pd.stable_entries
        
        for entry in all_entries_list:
            comp = entry.composition
            # 두 번째 원소의 composition fraction 계산
            if elements[1] in comp:
                fraction = comp.get_atomic_fraction(elements[1])
            else:
                fraction = 0.0
            
            # Formation energy per atom 계산 - PhaseDiagram의 reference state 기준으로 계산
            # pd.get_form_energy_per_atom은 PhaseDiagram이 생성될 때 설정된 reference state를 기준으로 계산
            # 이는 Materials Project와 동일한 방식입니다
            try:
                energy = pd.get_form_energy_per_atom(comp)
            except Exception as e:
                # get_form_energy_per_atom이 실패하면 entry를 사용하여 직접 계산
                # PhaseDiagram의 get_form_energy 메서드 사용 시도
                try:
                    energy = pd.get_form_energy(entry) / comp.num_atoms
                except:
                    # 최후의 수단: entry의 energy_per_atom 사용 (이것은 formation energy가 아님)
                    print(f"Warning: Could not calculate formation energy for {comp.reduced_formula}, using entry.energy_per_atom (not formation energy): {entry.energy_per_atom:.4f} eV/atom")
                    energy = entry.energy_per_atom
            
            formula = comp.reduced_formula
            
            # 안정적인 상인지 확인
            if entry in pd.stable_entries:
                stable_x.append(fraction)
                stable_y.append(energy)
                stable_labels.append(formula)
            else:
                # 불안정한 상 (above hull)
                # energy_above_hull 계산
                try:
                    e_above_hull = pd.get_e_above_hull(entry)
                except:
                    e_above_hull = 0.0  # 계산 실패 시 0으로 설정
                
                unstable_x.append(fraction)
                unstable_y.append(energy)
                unstable_labels.append(formula)
                unstable_e_above_hull.append(e_above_hull)
        
        # 안정적인 상 플롯
        if stable_x:
            ax.scatter(stable_x, stable_y, marker='o', 
                      edgecolor='black', facecolor='green', s=100, zorder=2, label='Stable')
            print(f"안정적인 상: {len(stable_x)}개")
        
        # 불안정한 상 플롯 (energy_above_hull에 따라 색상 변경)
        if unstable_x:
            # energy_above_hull 값을 0~1 범위로 정규화
            if unstable_e_above_hull:
                min_e_above_hull = min(unstable_e_above_hull)
                max_e_above_hull = max(unstable_e_above_hull)
                
                # 정규화 (0에 가까울수록 0, 1에 가까울수록 1)
                # 사용자가 "1에 가까울수록"이라고 했으므로, max_e_above_hull을 1로 매핑
                # 하지만 실제 값이 1보다 클 수 있으므로, max_e_above_hull을 기준으로 정규화
                if max_e_above_hull > min_e_above_hull:
                    normalized_values = [(e - min_e_above_hull) / (max_e_above_hull - min_e_above_hull) 
                                       for e in unstable_e_above_hull]
                else:
                    # 모든 값이 같으면 0으로 설정
                    normalized_values = [0.0] * len(unstable_e_above_hull)
                
                # YlOrRd 컬러맵 사용 (0에 가까울수록 Yellow, 1에 가까울수록 Red)
                from matplotlib.cm import get_cmap
                cmap = get_cmap('YlOrRd')
                colors = [cmap(val) for val in normalized_values]
            else:
                colors = 'red'  # energy_above_hull 값이 없으면 빨간색
            
            ax.scatter(unstable_x, unstable_y, marker='D', 
                      edgecolor='black', facecolor=colors, s=20, zorder=3, label='Above hull')
            print(f"불안정한 상 (above hull): {len(unstable_x)}개")
        
        # Convex hull 그리기 (안정적인 상들만)
        if len(stable_x) >= 2:
            # x 좌표로 정렬
            sorted_indices = sorted(range(len(stable_x)), key=lambda i: stable_x[i])
            sorted_x = [stable_x[i] for i in sorted_indices]
            sorted_y = [stable_y[i] for i in sorted_indices]
            
            # Reference state 인덱스 찾기 (composition fraction이 0 또는 1인 점)
            ref_state_indices = set()
            for i, x in enumerate(sorted_x):
                if abs(x) < 1e-6 or abs(x - 1.0) < 1e-6:
                    ref_state_indices.add(i)
            
            # Convex hull 계산
            from scipy.spatial import ConvexHull
            try:
                points = np.array([[x, y] for x, y in zip(sorted_x, sorted_y)])
                hull = ConvexHull(points)
                # Hull 경계 그리기 (reference state 두 개를 직접 연결하는 선 제외)
                for simplex in hull.simplices:
                    # Reference state 두 개를 직접 연결하는 선인지 확인
                    i1, i2 = simplex[0], simplex[1]
                    if i1 in ref_state_indices and i2 in ref_state_indices:
                        # Reference state 두 개를 직접 연결하는 선은 제외
                        continue
                    ax.plot(points[simplex, 0], points[simplex, 1], 
                           'k-', linewidth=1.5, alpha=0.7, zorder=1)
            except:
                # Convex hull 계산 실패 시 단순히 선으로 연결 (reference state 두 개 연결 제외)
                for i in range(len(sorted_x) - 1):
                    # Reference state 두 개를 직접 연결하는 선은 제외
                    if (i in ref_state_indices and (i + 1) in ref_state_indices):
                        continue
                    ax.plot([sorted_x[i], sorted_x[i+1]], [sorted_y[i], sorted_y[i+1]], 
                           'k-', linewidth=1.5, alpha=0.7, zorder=1)
        
        # Annotation 추가 (안정적인 상만)
        for i, (x, y, label) in enumerate(zip(stable_x, stable_y, stable_labels)):
            # 숫자를 아래첨자로 변환 (예: "Re2O3" → "Re$_2$O$_3$")
            # 화학식에서 숫자를 찾아서 LaTeX 아래첨자 형식으로 변환
            label_with_subscript = re.sub(r'(\d+)', r'$_{\1}$', label)
            
            # y=0이면 마커 위에, y<0이면 마커 오른쪽에 표시
            if abs(y) < 1e-6:  # y가 거의 0인 경우 (부동소수점 오차 고려)
                xytext = (0, 7)  # 마커 위쪽
                ha = 'center'  # 수평 중앙 정렬
                va = 'bottom'  # 수직 하단 정렬
            elif x > 0.5:
                xytext = (10, 0)  # 마커 오른쪽
                ha = 'left'  # 수평 중앙 정렬
                va = 'center'  # 수직 왼쪽 정렬
            elif x < 0.5:
                xytext = (0, 10)  # 기본값 (오른쪽 위)
                ha = 'right'  # 수평 중앙 정렬
                va = 'center'  # 수직 중앙 정렬
            
            ax.annotate(label_with_subscript, (x, y), xytext=xytext, 
                       textcoords='offset points', ha=ha, va=va,
                       fontsize=fontsize,
                       color='black',
                       fontweight='normal')
        
        # 축 설정
        ax.set_xlabel(f'Composition (Fraction {elements[1]})', 
                     fontsize=fontsize, fontweight='normal')
        ax.set_ylabel('Formation Energy (eV/atom)', 
                     fontsize=fontsize, fontweight='normal')
        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        
        # 그리드 추가
        ax.grid(True, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.5)
        ax.set_axisbelow(True)
        
        # 범위 설정
        if xlim is not None:
            ax.set_xlim(xlim)
        else:
            # 자동 범위 계산
            all_x = stable_x + unstable_x
            if all_x:
                x_margin = (max(all_x) - min(all_x)) * 0.1
                ax.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)
        
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            # 자동 범위 계산
            all_y = stable_y + unstable_y
            if all_y:
                y_margin = (max(all_y) - min(all_y)) * 0.1
                ax.set_ylim(min(all_y) - y_margin, max(all_y) + y_margin)
        
        print(f"X축 범위: {ax.get_xlim()[0]:.4f} ~ {ax.get_xlim()[1]:.4f}")
        print(f"Y축 범위: {ax.get_ylim()[0]:.4f} ~ {ax.get_ylim()[1]:.4f} eV/atom")
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Phase diagram이 {output_file}에 저장되었습니다.")
        
        if show_plot:
            plt.show()
            
    elif num_elements == 3:
        # 3원계: 삼각형 plot
        print("3원계 phase diagram을 그리는 중...")
        plotter = PDPlotter(pd, show_unstable=True, backend="matplotlib")
        plot = plotter.get_ternary_plot()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Phase diagram이 {output_file}에 저장되었습니다.")
        
        if show_plot:
            plt.show()
            
    else:
        print(f"{num_elements}원계는 시각화가 복잡합니다. 2원계 또는 3원계를 권장합니다.")
        # 4원계 이상은 projection을 사용해야 함
        plotter = PDPlotter(pd, show_unstable=True, backend="matplotlib")
        plot = plotter.get_plot()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Phase diagram이 {output_file}에 저장되었습니다.")
        
        if show_plot:
            plt.show()


def print_stable_phases(pd):
    """
    안정적인 상들을 출력합니다.
    
    Args:
        pd: PhaseDiagram 객체
    """
    print("\n=== 안정적인 상 (Stable Phases) ===")
    for i, entry in enumerate(pd.stable_entries, 1):
        formula = entry.composition.reduced_formula
        energy = entry.energy_per_atom
        print(f"{i:2d}. {formula:15s} - Energy: {energy:.4f} eV/atom")


def main():
    """
    메인 함수
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Materials Project Phase Diagram Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예제:
  # Re-O 2원계 phase diagram
  python Re-phase-diagram.py Re O
  
  # Re-Ru-O 3원계 phase diagram
  python Re-phase-diagram.py Re Ru O
  
  # 출력 파일 지정
  python Re-phase-diagram.py Re O -o re_o_phase_diagram.png
  
  # API 키 직접 지정
  python Re-phase-diagram.py Re O --api-key YOUR_API_KEY
  
  # 그림 사이즈 지정
  python Re-phase-diagram.py Re O --figsize 8 6
        """
    )
    
    parser.add_argument('elements', nargs='+', 
                       help='원소 기호들 (예: Re O 또는 Re Ru O)')
    parser.add_argument('-o', '--output', type=str, default=None,
                       help='출력 파일 경로 (예: phase_diagram.png)')
    parser.add_argument('--api-key', type=str, default=None,
                       help='Materials Project API 키 (없으면 MAPI_KEY 환경변수 사용)')
    parser.add_argument('--no-show', action='store_true',
                       help='플롯을 화면에 표시하지 않음')
    parser.add_argument('--list-stable', action='store_true',
                       help='안정적인 상 목록 출력')
    parser.add_argument('--xlim', type=float, nargs=2, metavar=('XMIN', 'XMAX'),
                       help='X축 범위 지정 (예: --xlim 0 1)')
    parser.add_argument('--ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'),
                       help='Y축 범위 지정 (예: --ylim -10 0)')
    parser.add_argument('--figsize', type=float, nargs=2, metavar=('WIDTH', 'HEIGHT'),
                       help='그림 사이즈 지정 인치 단위 (예: --figsize 8 6)')
    parser.add_argument('--stable-marker', type=str, default=None,
                       help='안정적인 상의 마커 스타일 (예: --stable-marker o 또는 s, ^, v, d 등)')
    parser.add_argument('--stable-color', type=str, default=None,
                       help='안정적인 상의 마커 색상 (예: --stable-color blue 또는 #0000FF)')
    parser.add_argument('--unstable-marker', type=str, default=None,
                       help='불안정한 상의 마커 스타일 (예: --unstable-marker s 또는 o, ^, v, d 등)')
    parser.add_argument('--unstable-color', type=str, default=None,
                       help='불안정한 상의 마커 색상 (예: --unstable-color red 또는 #FF0000)')
    parser.add_argument('--fontsize', type=int, default=14,
                       help='텍스트 크기 지정 (예: --fontsize 14)')
    
    args = parser.parse_args()
    
    # 원소 리스트
    elements = [e.capitalize() for e in args.elements]
    
    print(f"원소 시스템: {'-'.join(elements)}")
    print("=" * 50)
    
    # Phase diagram 가져오기
    pd, entries = get_phase_diagram(elements, api_key=args.api_key)
    
    # 안정적인 상 목록 출력
    if args.list_stable:
        print_stable_phases(pd)
    
    # Phase diagram 그리기
    output_file = args.output or f"{'_'.join(elements)}_phase_diagram.png"
    xlim = tuple(args.xlim) if args.xlim else None
    ylim = tuple(args.ylim) if args.ylim else None
    figsize = tuple(args.figsize) if args.figsize else None
    plot_phase_diagram(pd, elements, entries=entries, output_file=output_file, 
                      show_plot=not args.no_show, xlim=xlim, ylim=ylim, figsize=figsize, fontsize=args.fontsize)


if __name__ == "__main__":
    main()

