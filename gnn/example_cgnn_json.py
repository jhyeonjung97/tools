"""
Crystal Graph Neural Network (CGNN) - JSON 파일 사용 예제

CGNN의 핵심:
1. 각 노드(원자)에 대한 정보를 직접 넣을 수 있음
   - Bader charge, magnetic moment, coordination number 등
2. 노드 간 그래프 관계를 이용
   - 원자 간 거리 기반 엣지로 결정 구조 표현
"""

import torch
import torch.nn as nn
import torch.optim as optim
from gnn import CrystalGraphNeuralNetwork, json_to_graph, batch_json_to_graphs
import numpy as np
import json
from ase import Atoms
from ase.io import write
import os


def create_example_json_with_features():
    """각 원자별 정보를 포함한 예제 JSON 파일 생성"""
    print("=" * 60)
    print("예제 JSON 파일 생성 (각 원자별 정보 포함)")
    print("=" * 60)
    
    # 간단한 결정 구조 생성 (예: TiO2)
    atoms = Atoms('TiO2', 
                  positions=[[0, 0, 0], [2, 2, 2], [1, 1, 1]],
                  cell=[5, 5, 5],
                  pbc=[True, True, True])
    
    # 각 원자별 정보 추가
    num_atoms = len(atoms)
    
    # Bader charge (각 원자별)
    bader_charges = np.array([1.2, -0.6, -0.6])  # Ti, O, O
    
    # Magnetic moment (각 원자별)
    magmoms = np.array([0.0, 0.0, 0.0])
    
    # Coordination number (각 원자별)
    coordinations = np.array([6, 3, 3])  # Ti는 6배위, O는 3배위
    
    # ASE arrays에 저장
    atoms.arrays['bader_charge'] = bader_charges
    atoms.arrays['magmom'] = magmoms
    atoms.arrays['coordination'] = coordinations
    
    # JSON 파일로 저장
    output_file = 'example_structure.json'
    write(output_file, atoms)
    
    print(f"\n{output_file} 파일 생성 완료")
    print(f"\n각 원자별 정보:")
    print(f"  - Bader charge: {bader_charges}")
    print(f"  - Magnetic moment: {magmoms}")
    print(f"  - Coordination: {coordinations}")
    
    return output_file


def example_cgnn_with_custom_features():
    """커스텀 특징을 포함한 CGNN 사용 예제"""
    print("\n" + "=" * 60)
    print("커스텀 특징을 포함한 CGNN 사용 예제")
    print("=" * 60)
    
    json_file = "example_structure.json"
    
    if not os.path.exists(json_file):
        print(f"\n{json_file} 파일이 없습니다.")
        print("예제 파일을 생성합니다...")
        json_file = create_example_json_with_features()
    
    # JSON 파일 읽기 (각 원자별 정보 포함)
    print(f"\n{json_file} 파일 읽는 중...")
    data = json_to_graph(json_file, cutoff=5.0,
                        extra_feature_keys=['bader_charge', 'magmom', 'coordination'])
    
    print(f"\n그래프 정보:")
    print(f"  - 노드 개수: {data.x.shape[0]}")
    print(f"  - 노드 특징 차원: {data.x.shape[1]}")
    print(f"    * 원소 특징: {len(get_atomic_features(22, 'onehot'))}")  # Ti 원자 번호
    print(f"    * 추가 특징: Bader charge, Magmom, Coordination")
    print(f"  - 엣지 개수: {data.edge_index.shape[1] // 2}")
    
    # CGNN 모델 생성
    model = CrystalGraphNeuralNetwork(
        num_features=data.x.shape[1],  # 원소 특징 + 추가 특징들
        hidden_dim=64,
        num_outputs=1,  # 형성 에너지 등 예측
        num_layers=3,
        pooling='attention',  # Attention 기반 pooling
        use_edge_attr=True
    )
    
    print(f"\n모델 정보:")
    print(f"  - 입력 특징 차원: {data.x.shape[1]} (원소 + 각 원자별 정보)")
    print(f"  - 파라미터 개수: {sum(p.numel() for p in model.parameters())}")
    
    # 예측
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index, data.edge_attr)
        print(f"\n예측 결과:")
        print(f"  형성 에너지 (예측): {output.item():.4f} eV")
    
    return model, data


def example_batch_processing_with_features():
    """여러 JSON 파일을 일괄 처리하는 예제"""
    print("\n" + "=" * 60)
    print("여러 JSON 파일 일괄 처리 예제")
    print("=" * 60)
    
    # 현재 디렉토리의 JSON 파일 찾기
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    
    if len(json_files) == 0:
        print("\n현재 디렉토리에 JSON 파일이 없습니다.")
        print("\n사용법:")
        print("""
# 여러 JSON 파일을 일괄 처리
json_files = ['structure1.json', 'structure2.json', ...]
graphs = batch_json_to_graphs(
    json_files, 
    cutoff=5.0,
    extra_feature_keys=['bader_charge', 'magmom', 'coordination']
)

# 각 그래프의 특징 차원 확인
for i, graph in enumerate(graphs):
    print(f"Graph {i}: {graph.x.shape[1]} features")
        """)
        return
    
    print(f"\n{len(json_files)}개의 JSON 파일 발견")
    
    # 일괄 처리
    graphs = batch_json_to_graphs(
        json_files[:10],  # 처음 10개만 처리
        cutoff=5.0,
        extra_feature_keys=['bader_charge', 'magmom']
    )
    
    print(f"\n처리 완료: {len(graphs)}개 그래프 생성")
    
    if len(graphs) > 0:
        feature_dims = [g.x.shape[1] for g in graphs]
        print(f"\n특징 차원 통계:")
        print(f"  - 평균: {np.mean(feature_dims):.1f}")
        print(f"  - 최소: {min(feature_dims)}")
        print(f"  - 최대: {max(feature_dims)}")
        print(f"\n  각 그래프는 원소 특징 + 추가 특징들을 포함합니다")


def example_training_with_custom_features():
    """커스텀 특징을 포함한 학습 예제"""
    print("\n" + "=" * 60)
    print("커스텀 특징을 포함한 CGNN 학습 예제")
    print("=" * 60)
    
    print("\nCGNN 학습의 핵심:")
    print("  1. 각 원자별 정보를 노드 특징으로 사용")
    print("  2. 원자 간 거리 기반 엣지로 그래프 관계 표현")
    print("  3. 그래프 레벨 집계로 전체 결정 특성 예측")
    
    print("\n예제 코드:")
    print("""
from gnn import CrystalGraphNeuralNetwork, batch_json_to_graphs
import torch
import torch.nn as nn
import torch.optim as optim

# 1. 데이터 준비 (각 원자별 정보 포함)
json_files = ['structure1.json', 'structure2.json', ...]
graphs = batch_json_to_graphs(
    json_files, 
    cutoff=5.0,
    extra_feature_keys=['bader_charge', 'magmom', 'coordination']
)

# 2. 레이블 준비 (각 구조의 형성 에너지 등)
labels = torch.tensor([-5.2, -4.8, ...])  # 각 구조의 형성 에너지

# 3. CGNN 모델 생성
model = CrystalGraphNeuralNetwork(
    num_features=graphs[0].x.shape[1],  # 원소 특징 + 각 원자별 정보
    hidden_dim=64,
    num_outputs=1,  # 형성 에너지 예측
    num_layers=3,
    pooling='attention'
)

# 4. 학습 루프
optimizer = optim.Adam(model.parameters(), lr=0.001)
criterion = nn.MSELoss()

model.train()
for epoch in range(100):
    total_loss = 0
    for graph, label in zip(graphs, labels):
        optimizer.zero_grad()
        output = model(graph.x, graph.edge_index, graph.edge_attr)
        loss = criterion(output.squeeze(), label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    
    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1}: Loss = {total_loss/len(graphs):.4f}")
    """)


def get_atomic_features(atomic_number, feature_type='onehot'):
    """원자 특징 생성 함수 (예제용)"""
    if feature_type == 'onehot':
        common_elements = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                          19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
        if atomic_number in common_elements:
            features = np.zeros(len(common_elements))
            features[common_elements.index(atomic_number)] = 1.0
        else:
            features = np.zeros(len(common_elements))
            features[-1] = 1.0
        return features
    return np.array([atomic_number / 118.0])


if __name__ == "__main__":
    # 예제 JSON 파일 생성
    create_example_json_with_features()
    
    # 커스텀 특징을 포함한 CGNN 사용
    example_cgnn_with_custom_features()
    
    # 일괄 처리 예제
    example_batch_processing_with_features()
    
    # 학습 예제
    example_training_with_custom_features()
    
    print("\n" + "=" * 60)
    print("예제 실행 완료!")
    print("=" * 60)
    print("\nCGNN의 핵심:")
    print("  ✅ 각 노드(원자)에 대한 정보를 직접 넣을 수 있음")
    print("  ✅ 노드 간 그래프 관계를 이용하여 결정 구조 표현")
    print("  ✅ JSON 파일이 각 원자별 정보 추가에 가장 적합함")
