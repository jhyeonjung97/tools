#!/usr/bin/env python3
"""
Ru-Re-O Ternary Phase Diagram Generator

이 스크립트는 DFT 계산된 에너지를 사용하여 Ru-Re-O 삼원계 상평형도를 생성합니다.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
import pandas as pd
import json
from typing import Dict, List, Tuple, Optional
import argparse
import os

class TernaryPhaseDiagram:
    """
    Ru-Re-O 삼원계 상평형도 생성 클래스
    """
    
    def __init__(self, energy_data: Dict[str, float]):
        """
        Args:
            energy_data: 물질 이름을 키로 하고 형성에너지를 값으로 하는 딕셔너리
                        예: {'RuO2': -2.5, 'ReO3': -1.8, 'RuReO4': -3.2}
        """
        self.energy_data = energy_data
        self.elements = ['Ru', 'Re', 'O']
        self.formation_energies = {}
        
    def calculate_formation_energy(self, compound: str, composition: Dict[str, int], 
                                 total_energy: float) -> float:
        """
        형성에너지 계산
        
        Args:
            compound: 화합물 이름
            composition: 원소별 원자 수 {'Ru': 1, 'Re': 1, 'O': 4}
            total_energy: DFT 계산된 총 에너지 (eV)
            
        Returns:
            형성에너지 (eV/atom)
        """
        # 참조 원소들의 에너지 (eV/atom)
        reference_energies = {
            'Ru': -8.243244,   # Ru metal
            'Re': -12.42632685,   # Re metal  
            'O': -10.5835884822723/2 
        }
        
        # 참조 에너지 계산
        ref_energy = sum(composition[element] * reference_energies[element] 
                        for element in self.elements if element in composition)
        
        # 총 원자 수
        total_atoms = sum(composition.values())
        
        # 형성에너지 = (총에너지 - 참조에너지) / 총원자수
        formation_energy = (total_energy - ref_energy) / total_atoms
        
        return formation_energy
    
    def composition_to_ternary_coords(self, composition: Dict[str, int]) -> Tuple[float, float, float]:
        """
        조성을 삼원계 좌표로 변환
        
        Args:
            composition: 원소별 원자 수
            
        Returns:
            (x, y, z) 삼원계 좌표 (x + y + z = 1)
        """
        total_atoms = sum(composition.values())
        
        x = composition.get('Ru', 0) / total_atoms
        y = composition.get('Re', 0) / total_atoms  
        z = composition.get('O', 0) / total_atoms
        
        return x, y, z
    
    def ternary_to_cartesian(self, x: float, y: float, z: float) -> Tuple[float, float]:
        """
        삼원계 좌표를 카르테시안 좌표로 변환
        
        Args:
            x, y, z: 삼원계 좌표 (x + y + z = 1)
            
        Returns:
            (x_cart, y_cart) 카르테시안 좌표
        """
        # 삼각형 꼭짓점 좌표 (Re: 위쪽, O: 오른쪽 아래, Ru: 왼쪽 아래)
        A = (0, 0)  # Ru (왼쪽 아래)
        B = (0.5, np.sqrt(3)/2)  # Re (위쪽)  
        C = (1, 0)  # O (오른쪽 아래)
        
        # 변환 공식
        x_cart = x * A[0] + y * B[0] + z * C[0]
        y_cart = x * A[1] + y * B[1] + z * C[1]
        
        return x_cart, y_cart
    
    def create_ternary_diagram(self, compounds: List[Dict], 
                             show_convex_hull: bool = True,
                             show_formation_energies: bool = True,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        삼원계 상평형도 생성
        
        Args:
            compounds: 화합물 정보 리스트
                      [{'name': 'RuO2', 'composition': {'Ru': 1, 'O': 2}, 'energy': -10.5}, ...]
            show_convex_hull: Convex hull 표시 여부
            show_formation_energies: 형성에너지 표시 여부
            save_path: 저장 경로
            
        Returns:
            matplotlib Figure 객체
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # 삼각형 그리기 (Re: 위쪽, O: 오른쪽 아래, Ru: 왼쪽 아래)
        triangle_points = [(0, 0), (0.5, np.sqrt(3)/2), (1, 0), (0, 0)]
        triangle_x = [p[0] for p in triangle_points]
        triangle_y = [p[1] for p in triangle_points]
        ax.plot(triangle_x, triangle_y, 'k-', linewidth=2)
        
        # 꼭짓점 라벨
        ax.text(0, -0.02, 'Ru', fontsize=14, ha='center', va='top')
        ax.text(0.5, np.sqrt(3)/2 + 0.02, 'Re', fontsize=14, ha='center', va='bottom')  
        ax.text(1, -0.02, 'O', fontsize=14, ha='center', va='top')
        
        # 화합물들 플롯
        points = []
        labels = []
        energies = []
        
        for compound in compounds:
            name = compound['name']
            composition = compound['composition']
            energy = compound.get('energy', 0)
            
            # 삼원계 좌표로 변환
            x, y, z = self.composition_to_ternary_coords(composition)
            x_cart, y_cart = self.ternary_to_cartesian(x, y, z)
            
            points.append([x_cart, y_cart])
            labels.append(name)
            energies.append(energy)
            
            # # 점 플롯
            # ax.scatter(x_cart, y_cart, s=100, c='red', alpha=0.7, zorder=5)
            
            # 라벨 추가
            if name == 'Ru' or name == 'Re' or name == 'O2':
                continue
            elif name == 'RuO2' or name == 'RuO4':
                ax.annotate(name, (x_cart, y_cart), xytext=(0, -17), 
                        textcoords='offset points', fontsize=10, ha='center')
            else:
                ax.annotate(name, (x_cart, y_cart), xytext=(5, 5), 
                        textcoords='offset points', fontsize=10, ha='left')
        
        # # Convex hull 그리기
        # if show_convex_hull and len(points) >= 3:
        #     points_array = np.array(points)
        #     try:
        #         hull = ConvexHull(points_array)
        #         hull_points = points_array[hull.vertices]
        #         hull_points = np.vstack([hull_points, hull_points[0]])  # 닫기
                
        #         ax.plot(hull_points[:, 0], hull_points[:, 1], 'b--', alpha=0.7, linewidth=1)
                
        #         # Hull 내부 채우기
        #         polygon = Polygon(hull_points[:-1], alpha=0.2, facecolor='blue', edgecolor='blue')
        #         ax.add_patch(polygon)
        #     except:
        #         print("Convex hull을 계산할 수 없습니다.")
        
        # 형성에너지 컬러맵
        if show_formation_energies and len(energies) > 0:
            scatter = ax.scatter([p[0] for p in points], [p[1] for p in points], 
                               c=energies, s=100, cmap='viridis', alpha=0.8, zorder=5)
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
            cbar.set_label('Formation Energy (eV/atom)', fontsize=12)
        
        # 축 설정
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.0)
        ax.set_aspect('equal')
        ax.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"다이어그램이 {save_path}에 저장되었습니다.")
        
        return fig

def load_energy_data(file_path: str) -> Dict[str, float]:
    """
    에너지 데이터 파일 로드
    
    Args:
        file_path: JSON 또는 CSV 파일 경로
        
    Returns:
        에너지 데이터 딕셔너리
    """
    if file_path.endswith('.json'):
        with open(file_path, 'r') as f:
            return json.load(f)
    elif file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
        return dict(zip(df['compound'], df['energy']))
    else:
        raise ValueError("지원하지 않는 파일 형식입니다. JSON 또는 CSV 파일을 사용하세요.")

def main():
    """
    메인 함수 - 명령행 인터페이스
    """
    parser = argparse.ArgumentParser(description='Ru-Re-O 삼원계 상평형도 생성')
    parser.add_argument('--data', type=str, help='에너지 데이터 파일 경로 (JSON/CSV)')
    parser.add_argument('--output', type=str, default='ternary_phase_diagram.png', 
                       help='출력 파일 경로')
    parser.add_argument('--no-hull', action='store_true', help='Convex hull 표시 안함')
    parser.add_argument('--no-energy', action='store_true', help='형성에너지 컬러맵 표시 안함')
    
    args = parser.parse_args()
    
    # DFT 계산된 에너지 데이터 (eV 단위)
    # 여기에 실제 DFT 계산값을 입력하세요
    example_compounds = [
        {'name': 'RuO2', 'composition': {'Ru': 1, 'O': 2}, 'energy': -19.935983625},
        # {'name': 'RuO2', 'composition': {'Ru': 1, 'O': 2}, 'energy': -19.787085315},
        {'name': 'RuO4', 'composition': {'Ru': 1, 'O': 4}, 'energy': -30.11586173},
        {'name': 'ReO2', 'composition': {'Re': 1, 'O': 2}, 'energy': -26.6135512425},
        {'name': 'ReO3', 'composition': {'Re': 1, 'O': 3}, 'energy': -33.65095848},
        {'name': 'Re2O7', 'composition': {'Re': 2, 'O': 7}, 'energy': -72.63795125},
        {'name': 'Re-RuO2', 'composition': {'Re': 1, 'Ru': 7, 'O': 16}, 'energy': -166.111256},
        {'name': 'Ru', 'composition': {'Ru': 1}, 'energy': -8.243244},
        {'name': 'Re', 'composition': {'Re': 1}, 'energy': -12.42632685},
        {'name': 'O2', 'composition': {'O': 2}, 'energy': -10.5835884822723},
    ]

    # 데이터 로드
    if args.data and os.path.exists(args.data):
        energy_data = load_energy_data(args.data)
        # 에너지 데이터를 화합물 정보로 변환하는 로직 추가 필요
        compounds = example_compounds  # 임시로 예시 데이터 사용
    else:
        compounds = example_compounds
    
    # 삼원계 상평형도 생성
    diagram = TernaryPhaseDiagram({})
    fig = diagram.create_ternary_diagram(
        compounds=compounds,
        show_convex_hull=not args.no_hull,
        show_formation_energies=not args.no_energy,
        save_path=args.output
    )
    
    plt.show()

if __name__ == "__main__":
    main()
