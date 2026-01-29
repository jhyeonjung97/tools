"""
JSON 파일을 GNN 입력으로 사용하는 예제
Bader charge 등 추가 정보를 포함한 그래프 생성
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from gnn import json_to_graph, batch_json_to_graphs
import os
import numpy as np


class CrystalGCNWithExtraFeatures(nn.Module):
    """
    추가 특징(Bader charge 등)을 포함한 결정 구조용 GCN 모델
    """
    def __init__(self, num_features, hidden_dim, num_classes, num_layers=3, dropout=0.5):
        super(CrystalGCNWithExtraFeatures, self).__init__()
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(num_features, hidden_dim))
        
        for _ in range(num_layers - 2):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
        
        self.convs.append(GCNConv(hidden_dim, num_classes))
        self.dropout = dropout
    
    def forward(self, x, edge_index):
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        x = self.convs[-1](x, edge_index)
        return F.log_softmax(x, dim=1)


def example_single_json():
    """단일 JSON 파일 처리 예제"""
    print("=" * 60)
    print("단일 JSON 파일 처리 예제")
    print("=" * 60)
    
    json_file = "example.json"
    
    if not os.path.exists(json_file):
        print(f"\n{json_file} 파일이 없습니다.")
        print("\n사용법:")
        print("  from gnn import json_to_graph")
        print("  ")
        print("  # 기본 사용")
        print("  data = json_to_graph('structure.json', cutoff=5.0)")
        print("  ")
        print("  # Bader charge 등 추가 특징 포함")
        print("  data = json_to_graph('structure.json', cutoff=5.0,")
        print("                       extra_feature_keys=['bader_charge', 'magmom'])")
        print("  ")
        print("  print(f'노드 개수: {data.x.shape[0]}')")
        print("  print(f'노드 특징 차원: {data.x.shape[1]}')")
        return
    
    # JSON 파일을 그래프로 변환
    print(f"\n{json_file} 파일 읽는 중...")
    
    # Bader charge 등 추가 특징 포함
    data = json_to_graph(json_file, cutoff=5.0, feature_type='onehot', 
                        extra_feature_keys=['bader_charge', 'magmom', 'forces'])
    
    print(f"\n그래프 정보:")
    print(f"  - 노드 개수: {data.x.shape[0]}")
    print(f"  - 노드 특징 차원: {data.x.shape[1]}")
    print(f"  - 엣지 개수: {data.edge_index.shape[1] // 2}")
    if hasattr(data, 'edge_attr'):
        print(f"  - 엣지 속성 차원: {data.edge_attr.shape[1]}")
    
    # 원자 종류 정보
    if hasattr(data, 'atomic_numbers'):
        unique_elements = torch.unique(data.atomic_numbers).tolist()
        print(f"  - 원소 종류: {unique_elements}")
    
    # 모델 생성 및 예측
    num_classes = 2  # 예제용
    model = CrystalGCNWithExtraFeatures(data.x.shape[1], hidden_dim=64, num_classes=num_classes)
    
    print(f"\n모델 정보:")
    print(f"  - 입력 특징 차원: {data.x.shape[1]}")
    print(f"  - 파라미터 개수: {sum(p.numel() for p in model.parameters())}")
    
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index)
        print(f"\n예측 결과:")
        print(f"  출력 형태: {output.shape}")
        print(f"  예측 클래스: {output.argmax(dim=1).tolist()}")
    
    return data, model


def example_batch_json():
    """여러 JSON 파일 일괄 처리 예제"""
    print("\n" + "=" * 60)
    print("여러 JSON 파일 일괄 처리 예제")
    print("=" * 60)
    
    # 현재 디렉토리의 모든 JSON 파일 찾기
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    
    if len(json_files) == 0:
        print("\n현재 디렉토리에 JSON 파일이 없습니다.")
        print("\n사용법:")
        print("  json_files = ['file1.json', 'file2.json', 'file3.json']")
        print("  graphs = batch_json_to_graphs(json_files, cutoff=5.0,")
        print("                                extra_feature_keys=['bader_charge'])")
        print("  print(f'처리된 그래프 개수: {len(graphs)}')")
        return
    
    print(f"\n{len(json_files)}개의 JSON 파일 발견:")
    for f in json_files[:5]:  # 처음 5개만 표시
        print(f"  - {f}")
    if len(json_files) > 5:
        print(f"  ... 외 {len(json_files) - 5}개")
    
    # 일괄 처리 (Bader charge 포함)
    print("\n그래프로 변환 중 (Bader charge 등 추가 특징 포함)...")
    graphs = batch_json_to_graphs(json_files[:10], cutoff=5.0, 
                                  extra_feature_keys=['bader_charge', 'magmom'])
    
    print(f"\n처리 완료: {len(graphs)}개 그래프 생성")
    
    # 통계 정보
    if len(graphs) > 0:
        num_nodes_list = [g.x.shape[0] for g in graphs]
        num_edges_list = [g.edge_index.shape[1] // 2 for g in graphs]
        feature_dims = [g.x.shape[1] for g in graphs]
        
        print(f"\n그래프 통계:")
        print(f"  - 평균 노드 개수: {sum(num_nodes_list) / len(num_nodes_list):.1f}")
        print(f"  - 평균 엣지 개수: {sum(num_edges_list) / len(num_edges_list):.1f}")
        print(f"  - 평균 특징 차원: {sum(feature_dims) / len(feature_dims):.1f}")
        print(f"  - 최소 노드 개수: {min(num_nodes_list)}")
        print(f"  - 최대 노드 개수: {max(num_nodes_list)}")
        print(f"  - 최소 특징 차원: {min(feature_dims)}")
        print(f"  - 최대 특징 차원: {max(feature_dims)}")
    
    return graphs


def example_create_json_with_extra_features():
    """추가 특징을 포함한 JSON 파일 생성 예제"""
    print("\n" + "=" * 60)
    print("추가 특징을 포함한 JSON 파일 생성 예제")
    print("=" * 60)
    
    if not os.path.exists("example_structure.json"):
        print("\n예제: ASE를 사용하여 추가 특징을 포함한 JSON 파일 생성")
        print("""
from ase import Atoms
from ase.io import write
import numpy as np

# 구조 생성
atoms = Atoms('H2O', positions=[[0, 0, 0], [0, 1, 0], [0, 0.5, 0.7]])

# 추가 특징 추가 (Bader charge 예시)
bader_charges = np.array([-0.5, 0.25, 0.25])  # 각 원자별 Bader charge
atoms.arrays['bader_charge'] = bader_charges

# 자기 모멘트 추가
magmoms = np.array([0.0, 0.0, 0.0])
atoms.arrays['magmom'] = magmoms

# JSON 파일로 저장
write('structure_with_features.json', atoms)

# 또는 커스텀 JSON 형식으로 저장
import json
custom_data = {
    'atoms': atoms,  # ASE atoms 객체
    'extra': {
        'bader_charge': bader_charges.tolist(),
        'magmom': magmoms.tolist(),
        'custom_feature': [1.0, 2.0, 3.0]  # 다른 커스텀 특징
    }
}
# JSON으로 저장하려면 atoms를 먼저 변환해야 함
        """)
        return
    
    print("\n예제 JSON 파일이 생성되었습니다.")


def example_training_with_json():
    """JSON 파일로 학습하는 예제"""
    print("\n" + "=" * 60)
    print("JSON 파일로 학습하는 예제")
    print("=" * 60)
    
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    
    if len(json_files) < 2:
        print("\n학습을 위해서는 최소 2개의 JSON 파일이 필요합니다.")
        print("\n예제 코드:")
        print("""
# 1. 데이터 준비
json_files = ['structure1.json', 'structure2.json', ...]
graphs = batch_json_to_graphs(json_files, cutoff=5.0,
                             extra_feature_keys=['bader_charge', 'magmom'])

# 2. 레이블 준비 (예: 밴드갭, 형성 에너지 등)
labels = torch.tensor([0.5, 1.2, ...])  # 예제 레이블

# 3. 모델 생성
model = CrystalGCNWithExtraFeatures(
    num_features=graphs[0].x.shape[1],  # 원소 특징 + 추가 특징
    hidden_dim=64, 
    num_classes=1  # 회귀 작업의 경우
)

# 4. 학습 루프
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
criterion = nn.MSELoss()

for epoch in range(100):
    total_loss = 0
    for graph, label in zip(graphs, labels):
        optimizer.zero_grad()
        output = model(graph.x, graph.edge_index)
        loss = criterion(output.mean(), label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    print(f"Epoch {epoch+1}: Loss = {total_loss/len(graphs):.4f}")
        """)
        return
    
    print(f"\n{len(json_files)}개의 JSON 파일로 학습 예제를 실행합니다.")
    print("(실제 학습을 위해서는 레이블 데이터가 필요합니다)")


if __name__ == "__main__":
    # 단일 JSON 파일 예제
    example_single_json()
    
    # 일괄 처리 예제
    example_batch_json()
    
    # JSON 파일 생성 예제
    example_create_json_with_extra_features()
    
    # 학습 예제
    example_training_with_json()
    
    print("\n" + "=" * 60)
    print("예제 실행 완료!")
    print("=" * 60)
