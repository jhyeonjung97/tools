"""
CIF 파일을 GNN 입력으로 사용하는 예제
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from gnn import cif_to_graph, batch_cif_to_graphs
import os


class CrystalGCN(nn.Module):
    """
    결정 구조를 위한 GCN 모델
    """
    def __init__(self, num_features, hidden_dim, num_classes, num_layers=3, dropout=0.5):
        super(CrystalGCN, self).__init__()
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


def example_single_cif():
    """단일 CIF 파일 처리 예제"""
    print("=" * 60)
    print("단일 CIF 파일 처리 예제")
    print("=" * 60)
    
    # CIF 파일 경로 (예제)
    cif_file = "example.cif"
    
    if not os.path.exists(cif_file):
        print(f"\n{cif_file} 파일이 없습니다.")
        print("사용법:")
        print("  data = cif_to_graph('your_structure.cif', cutoff=5.0)")
        print("  print(f'노드 개수: {data.x.shape[0]}')")
        print("  print(f'엣지 개수: {data.edge_index.shape[1] // 2}')")
        print("  print(f'노드 특징 차원: {data.x.shape[1]}')")
        return
    
    # CIF 파일을 그래프로 변환
    print(f"\n{cif_file} 파일 읽는 중...")
    data = cif_to_graph(cif_file, cutoff=5.0, feature_type='onehot', use_periodic=True)
    
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
    model = CrystalGCN(data.x.shape[1], hidden_dim=64, num_classes=num_classes)
    
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index)
        print(f"\n예측 결과:")
        print(f"  출력 형태: {output.shape}")
        print(f"  예측 클래스: {output.argmax(dim=1).tolist()}")
    
    return data, model


def example_batch_cif():
    """여러 CIF 파일 일괄 처리 예제"""
    print("\n" + "=" * 60)
    print("여러 CIF 파일 일괄 처리 예제")
    print("=" * 60)
    
    # 현재 디렉토리의 모든 CIF 파일 찾기
    cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
    
    if len(cif_files) == 0:
        print("\n현재 디렉토리에 CIF 파일이 없습니다.")
        print("사용법:")
        print("  cif_files = ['file1.cif', 'file2.cif', 'file3.cif']")
        print("  graphs = batch_cif_to_graphs(cif_files, cutoff=5.0)")
        print("  print(f'처리된 그래프 개수: {len(graphs)}')")
        return
    
    print(f"\n{len(cif_files)}개의 CIF 파일 발견:")
    for f in cif_files[:5]:  # 처음 5개만 표시
        print(f"  - {f}")
    if len(cif_files) > 5:
        print(f"  ... 외 {len(cif_files) - 5}개")
    
    # 일괄 처리
    print("\n그래프로 변환 중...")
    graphs = batch_cif_to_graphs(cif_files[:10], cutoff=5.0)  # 처음 10개만 처리
    
    print(f"\n처리 완료: {len(graphs)}개 그래프 생성")
    
    # 통계 정보
    if len(graphs) > 0:
        num_nodes_list = [g.x.shape[0] for g in graphs]
        num_edges_list = [g.edge_index.shape[1] // 2 for g in graphs]
        
        print(f"\n그래프 통계:")
        print(f"  - 평균 노드 개수: {sum(num_nodes_list) / len(num_nodes_list):.1f}")
        print(f"  - 평균 엣지 개수: {sum(num_edges_list) / len(num_edges_list):.1f}")
        print(f"  - 최소 노드 개수: {min(num_nodes_list)}")
        print(f"  - 최대 노드 개수: {max(num_nodes_list)}")
    
    return graphs


def example_training_with_cif():
    """CIF 파일로 학습하는 예제"""
    print("\n" + "=" * 60)
    print("CIF 파일로 학습하는 예제")
    print("=" * 60)
    
    # 예제: 여러 CIF 파일을 읽어서 학습
    cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
    
    if len(cif_files) < 2:
        print("\n학습을 위해서는 최소 2개의 CIF 파일이 필요합니다.")
        print("\n예제 코드:")
        print("""
# 1. 데이터 준비
cif_files = ['structure1.cif', 'structure2.cif', ...]
graphs = batch_cif_to_graphs(cif_files, cutoff=5.0)

# 2. 레이블 준비 (예: 밴드갭, 형성 에너지 등)
labels = torch.tensor([0.5, 1.2, ...])  # 예제 레이블

# 3. 모델 생성
model = CrystalGCN(num_features=graphs[0].x.shape[1], 
                   hidden_dim=64, 
                   num_classes=1)  # 회귀 작업의 경우

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
    
    print(f"\n{len(cif_files)}개의 CIF 파일로 학습 예제를 실행합니다.")
    print("(실제 학습을 위해서는 레이블 데이터가 필요합니다)")


if __name__ == "__main__":
    # 단일 CIF 파일 예제
    example_single_cif()
    
    # 일괄 처리 예제
    example_batch_cif()
    
    # 학습 예제
    example_training_with_cif()
    
    print("\n" + "=" * 60)
    print("예제 실행 완료!")
    print("=" * 60)
