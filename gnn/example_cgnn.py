"""
Crystal Graph Neural Network (CGNN) 사용 예제

CGNN은 결정 구조의 특성을 예측하기 위한 특화된 GNN 모델입니다.
Materials Project 등에서 널리 사용되는 아키텍처입니다.
"""

import torch
import torch.nn as nn
import torch.optim as optim
from gnn import CrystalGraphNeuralNetwork, cif_to_graph, json_to_graph
import os


def example_cgnn_regression():
    """CGNN을 사용한 회귀 작업 예제 (형성 에너지 예측 등)"""
    print("=" * 60)
    print("CGNN 회귀 작업 예제")
    print("=" * 60)
    
    # 예제: CIF 파일이 있는 경우
    cif_file = "example.cif"
    
    if not os.path.exists(cif_file):
        print(f"\n{cif_file} 파일이 없습니다.")
        print("\n사용법:")
        print("  from gnn import CrystalGraphNeuralNetwork, cif_to_graph")
        print("  ")
        print("  # 데이터 준비")
        print("  data = cif_to_graph('structure.cif', cutoff=5.0, add_edge_attr=True)")
        print("  ")
        print("  # 모델 생성")
        print("  model = CrystalGraphNeuralNetwork(")
        print("      num_features=data.x.shape[1],")
        print("      hidden_dim=64,")
        print("      num_outputs=1  # 회귀 작업")
        print("  )")
        print("  ")
        print("  # 예측")
        print("  output = model(data.x, data.edge_index, data.edge_attr)")
        return
    
    # 데이터 준비
    print(f"\n{cif_file} 파일 읽는 중...")
    data = cif_to_graph(cif_file, cutoff=5.0, add_edge_attr=True)
    
    print(f"\n그래프 정보:")
    print(f"  - 노드 개수: {data.x.shape[0]}")
    print(f"  - 노드 특징 차원: {data.x.shape[1]}")
    print(f"  - 엣지 개수: {data.edge_index.shape[1] // 2}")
    
    # 모델 생성
    model = CrystalGraphNeuralNetwork(
        num_features=data.x.shape[1],
        hidden_dim=64,
        num_outputs=1,  # 회귀 작업 (형성 에너지 등)
        num_layers=3,
        pooling='mean',
        use_edge_attr=True
    )
    
    print(f"\n모델 정보:")
    print(f"  - 입력 특징 차원: {data.x.shape[1]}")
    print(f"  - 은닉 차원: 64")
    print(f"  - 출력 차원: 1")
    print(f"  - 파라미터 개수: {sum(p.numel() for p in model.parameters())}")
    
    # 예측
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index, data.edge_attr)
        print(f"\n예측 결과:")
        print(f"  예측값: {output.item():.4f}")
    
    return model, data


def example_cgnn_training():
    """CGNN 학습 예제"""
    print("\n" + "=" * 60)
    print("CGNN 학습 예제")
    print("=" * 60)
    
    # 예제: 여러 구조 파일이 있는 경우
    cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    
    if len(cif_files) < 2 and len(json_files) < 2:
        print("\n학습을 위해서는 최소 2개의 구조 파일이 필요합니다.")
        print("\n예제 코드:")
        print("""
from gnn import CrystalGraphNeuralNetwork, batch_cif_to_graphs
import torch
import torch.nn as nn
import torch.optim as optim

# 1. 데이터 준비
cif_files = ['structure1.cif', 'structure2.cif', ...]
graphs = batch_cif_to_graphs(cif_files, cutoff=5.0)

# 2. 레이블 준비 (예: 형성 에너지)
labels = torch.tensor([-5.2, -4.8, ...])  # 각 구조의 형성 에너지

# 3. 모델 생성
model = CrystalGraphNeuralNetwork(
    num_features=graphs[0].x.shape[1],
    hidden_dim=64,
    num_outputs=1,  # 회귀 작업
    num_layers=3,
    pooling='mean'
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
        return
    
    print(f"\n{len(cif_files)}개의 CIF 파일과 {len(json_files)}개의 JSON 파일 발견")
    print("(실제 학습을 위해서는 레이블 데이터가 필요합니다)")


def example_cgnn_with_json():
    """JSON 파일을 사용한 CGNN 예제"""
    print("\n" + "=" * 60)
    print("JSON 파일을 사용한 CGNN 예제")
    print("=" * 60)
    
    json_file = "example.json"
    
    if not os.path.exists(json_file):
        print(f"\n{json_file} 파일이 없습니다.")
        print("\n사용법:")
        print("  from gnn import CrystalGraphNeuralNetwork, json_to_graph")
        print("  ")
        print("  # JSON 파일에서 그래프 생성 (Bader charge 포함)")
        print("  data = json_to_graph('structure.json', cutoff=5.0,")
        print("                       extra_feature_keys=['bader_charge'])")
        print("  ")
        print("  # CGNN 모델 생성")
        print("  model = CrystalGraphNeuralNetwork(")
        print("      num_features=data.x.shape[1],  # 원소 특징 + Bader charge")
        print("      hidden_dim=64,")
        print("      num_outputs=1")
        print("  )")
        return
    
    # JSON 파일 읽기
    print(f"\n{json_file} 파일 읽는 중...")
    data = json_to_graph(json_file, cutoff=5.0, 
                        extra_feature_keys=['bader_charge', 'magmom'])
    
    print(f"\n그래프 정보:")
    print(f"  - 노드 개수: {data.x.shape[0]}")
    print(f"  - 노드 특징 차원: {data.x.shape[1]}")  # 원소 특징 + 추가 특징
    print(f"  - 엣지 개수: {data.edge_index.shape[1] // 2}")
    
    # 모델 생성
    model = CrystalGraphNeuralNetwork(
        num_features=data.x.shape[1],
        hidden_dim=64,
        num_outputs=1,
        num_layers=3,
        pooling='attention',  # Attention 기반 pooling 사용
        use_edge_attr=True
    )
    
    # 예측
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index, data.edge_attr)
        print(f"\n예측 결과:")
        print(f"  예측값: {output.item():.4f}")
    
    return model, data


def example_different_pooling():
    """다양한 pooling 방법 비교 예제"""
    print("\n" + "=" * 60)
    print("다양한 Pooling 방법 비교")
    print("=" * 60)
    
    print("\nCGNN은 다양한 그래프 집계 방법을 지원합니다:")
    print("  - 'mean': 평균 집계 (가장 일반적)")
    print("  - 'sum': 합산 집계")
    print("  - 'max': 최대값 집계")
    print("  - 'attention': Attention 기반 가중 평균 (가장 정교함)")
    
    print("\n사용 예:")
    print("  # Mean pooling")
    print("  model_mean = CrystalGraphNeuralNetwork(..., pooling='mean')")
    print("  ")
    print("  # Attention pooling")
    print("  model_att = CrystalGraphNeuralNetwork(..., pooling='attention')")
    
    print("\n일반적으로 attention pooling이 더 좋은 성능을 보이지만,")
    print("계산 비용이 더 높습니다.")


if __name__ == "__main__":
    # 회귀 작업 예제
    example_cgnn_regression()
    
    # 학습 예제
    example_cgnn_training()
    
    # JSON 파일 사용 예제
    example_cgnn_with_json()
    
    # Pooling 방법 비교
    example_different_pooling()
    
    print("\n" + "=" * 60)
    print("예제 실행 완료!")
    print("=" * 60)
