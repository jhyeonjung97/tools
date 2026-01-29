"""
Graph Neural Network (GNN) 구현 예제

이 모듈은 기본적인 Graph Convolutional Network (GCN) 레이어와 
Crystal Graph Neural Network (CGNN)을 구현합니다.

Crystal Graph Neural Network는 결정 구조의 특성을 예측하기 위한
특화된 GNN 모델로, Materials Project 등에서 널리 사용됩니다.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.datasets import Planetoid
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List
from collections import Counter

# CIF 파일 읽기를 위한 라이브러리
try:
    from ase.io import read as ase_read
    from ase import Atoms
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    print("Warning: ASE가 설치되지 않았습니다. CIF/JSON 파일 읽기 기능을 사용하려면 설치하세요: pip install ase")

try:
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core import Structure
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False

# JSON 파일 읽기를 위한 라이브러리
import json


class BasicGCNLayer(nn.Module):
    """
    기본적인 Graph Convolutional Network 레이어 구현
    
    수식: H^(l+1) = σ(D^(-1/2) A D^(-1/2) H^(l) W^(l))
    여기서:
    - A: 인접 행렬
    - D: 차수 행렬 (degree matrix)
    - H^(l): l번째 레이어의 노드 특징
    - W^(l): 학습 가능한 가중치 행렬
    """
    
    def __init__(self, in_features, out_features, bias=True):
        super(BasicGCNLayer, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()
    
    def reset_parameters(self):
        """가중치 초기화"""
        nn.init.xavier_uniform_(self.weight)
        if self.bias is not None:
            nn.init.zeros_(self.bias)
    
    def forward(self, x, adj):
        """
        Forward pass
        
        Args:
            x: 노드 특징 행렬 [N, in_features]
            adj: 정규화된 인접 행렬 [N, N]
        
        Returns:
            출력 특징 행렬 [N, out_features]
        """
        # 특징 변환: XW
        support = torch.mm(x, self.weight)
        # 그래프 컨볼루션: AXW
        output = torch.mm(adj, support)
        
        if self.bias is not None:
            output += self.bias
        
        return output


class SimpleGCN(nn.Module):
    """
    간단한 2-layer GCN 모델
    """
    
    def __init__(self, num_features, hidden_dim, num_classes, dropout=0.5):
        super(SimpleGCN, self).__init__()
        self.gcn1 = BasicGCNLayer(num_features, hidden_dim)
        self.gcn2 = BasicGCNLayer(hidden_dim, num_classes)
        self.dropout = dropout
    
    def forward(self, x, adj):
        """
        Forward pass
        
        Args:
            x: 노드 특징 [N, num_features]
            adj: 정규화된 인접 행렬 [N, N]
        
        Returns:
            로그 확률 [N, num_classes]
        """
        # 첫 번째 GCN 레이어
        x = self.gcn1(x, adj)
        x = F.relu(x)
        x = F.dropout(x, self.dropout, training=self.training)
        
        # 두 번째 GCN 레이어
        x = self.gcn2(x, adj)
        
        return F.log_softmax(x, dim=1)


class CrystalGraphConv(nn.Module):
    """
    Crystal Graph Convolution 레이어
    
    결정 구조에 특화된 메시지 전달 방식:
    - 엣지 속성(거리)을 사용한 가중 메시지 전달
    - 이웃 노드의 특징을 거리 기반으로 집계
    
    수식: h_i^(l+1) = σ(Σ_{j∈N(i)} W^(l) * h_j^(l) * f(d_ij) + b^(l))
    여기서 f(d_ij)는 거리 기반 가중치 함수
    """
    
    def __init__(self, in_features, out_features, edge_dim=1, bias=True):
        """
        Args:
            in_features: 입력 특징 차원
            out_features: 출력 특징 차원
            edge_dim: 엣지 특징 차원 (거리 등)
            bias: bias 사용 여부
        """
        super(CrystalGraphConv, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.edge_dim = edge_dim
        
        # 노드 특징 변환 가중치
        self.weight = nn.Parameter(torch.FloatTensor(in_features, out_features))
        
        # 엣지 특징(거리)을 사용한 가중치 네트워크
        if edge_dim > 0:
            self.edge_weight_net = nn.Sequential(
                nn.Linear(edge_dim, out_features),
                nn.Sigmoid()  # 거리 기반 가중치를 [0, 1] 범위로
            )
        
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        
        self.reset_parameters()
    
    def reset_parameters(self):
        """가중치 초기화"""
        nn.init.xavier_uniform_(self.weight)
        if self.bias is not None:
            nn.init.zeros_(self.bias)
    
    def forward(self, x, edge_index, edge_attr=None):
        """
        Forward pass
        
        Args:
            x: 노드 특징 [N, in_features]
            edge_index: 엣지 인덱스 [2, E]
            edge_attr: 엣지 속성 [E, edge_dim] (거리 등)
        
        Returns:
            출력 특징 [N, out_features]
        """
        # 노드 특징 변환
        x_transformed = torch.mm(x, self.weight)  # [N, out_features]
        
        # 메시지 전달을 위한 출력 텐서 초기화
        out = torch.zeros(x.size(0), self.out_features, device=x.device, dtype=x.dtype)
        
        # 엣지별 메시지 집계
        if edge_attr is not None and self.edge_dim > 0:
            # 엣지 가중치 계산 (거리 기반)
            edge_weights = self.edge_weight_net(edge_attr)  # [E, out_features]
            
            # 각 엣지에 대해 메시지 전달
            for k in range(edge_index.size(1)):
                i = edge_index[0, k]  # source node
                j = edge_index[1, k]  # target node
                # 메시지: 이웃 노드의 변환된 특징 * 거리 기반 가중치
                out[j] += x_transformed[i] * edge_weights[k]
        else:
            # 엣지 속성이 없는 경우 단순 합산
            for k in range(edge_index.size(1)):
                i = edge_index[0, k]
                j = edge_index[1, k]
                out[j] += x_transformed[i]
        
        # Bias 추가
        if self.bias is not None:
            out += self.bias
        
        return out


class CrystalGraphNeuralNetwork(nn.Module):
    """
    Crystal Graph Neural Network (CGNN)
    
    결정 구조의 특성을 예측하기 위한 특화된 GNN 모델
    - Crystal Graph Convolution 레이어 사용
    - 전체 그래프에 대한 집계 (global pooling)로 결정 특성 예측
    - 회귀 및 분류 작업 모두 지원
    """
    
    def __init__(self, 
                 num_features, 
                 hidden_dim, 
                 num_outputs=1,
                 num_layers=3,
                 dropout=0.5,
                 pooling='mean',
                 use_edge_attr=True):
        """
        Args:
            num_features: 입력 노드 특징 차원
            hidden_dim: 은닉층 차원
            num_outputs: 출력 차원 (회귀: 1, 분류: 클래스 수)
            num_layers: 레이어 개수
            dropout: Dropout 비율
            pooling: 그래프 집계 방법 ('mean', 'sum', 'max', 'attention')
            use_edge_attr: 엣지 속성 사용 여부
        """
        super(CrystalGraphNeuralNetwork, self).__init__()
        self.num_layers = num_layers
        self.dropout = dropout
        self.pooling = pooling
        self.use_edge_attr = use_edge_attr
        
        # Crystal Graph Convolution 레이어들
        self.convs = nn.ModuleList()
        
        # 첫 번째 레이어
        self.convs.append(CrystalGraphConv(num_features, hidden_dim, 
                                          edge_dim=1 if use_edge_attr else 0))
        
        # 중간 레이어들
        for _ in range(num_layers - 2):
            self.convs.append(CrystalGraphConv(hidden_dim, hidden_dim,
                                              edge_dim=1 if use_edge_attr else 0))
        
        # 마지막 레이어
        self.convs.append(CrystalGraphConv(hidden_dim, hidden_dim,
                                          edge_dim=1 if use_edge_attr else 0))
        
        # Attention 기반 pooling을 위한 네트워크
        if pooling == 'attention':
            self.attention = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.Tanh(),
                nn.Linear(hidden_dim, 1)
            )
        
        # 출력 레이어
        self.output_layer = nn.Linear(hidden_dim, num_outputs)
    
    def forward(self, x, edge_index, edge_attr=None, batch=None):
        """
        Forward pass
        
        Args:
            x: 노드 특징 [N, num_features]
            edge_index: 엣지 인덱스 [2, E]
            edge_attr: 엣지 속성 [E, 1] (거리 등)
            batch: 배치 인덱스 [N] (여러 그래프를 배치로 처리할 때)
        
        Returns:
            출력 [batch_size, num_outputs] 또는 [num_outputs]
        """
        # 엣지 속성 처리
        if edge_attr is None and self.use_edge_attr:
            # 엣지 속성이 없는 경우 거리 계산 (간단한 예제)
            # 실제로는 미리 계산된 거리를 사용하는 것이 좋음
            edge_attr = torch.ones(edge_index.size(1), 1, device=x.device)
        
        # Crystal Graph Convolution 레이어들 통과
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_attr if self.use_edge_attr else None)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        # 마지막 레이어
        x = self.convs[-1](x, edge_index, edge_attr if self.use_edge_attr else None)
        x = F.relu(x)
        
        # 그래프 레벨 집계 (Global Pooling)
        if batch is None:
            # 단일 그래프인 경우
            graph_embedding = self._pool(x)
        else:
            # 배치 처리
            graph_embedding = self._pool_batch(x, batch)
        
        # 출력
        output = self.output_layer(graph_embedding)
        
        return output
    
    def _pool(self, x):
        """단일 그래프에 대한 집계"""
        if self.pooling == 'mean':
            return x.mean(dim=0, keepdim=True)
        elif self.pooling == 'sum':
            return x.sum(dim=0, keepdim=True)
        elif self.pooling == 'max':
            return x.max(dim=0, keepdim=True)[0]
        elif self.pooling == 'attention':
            # Attention 기반 가중 평균
            att_weights = self.attention(x)  # [N, 1]
            att_weights = F.softmax(att_weights, dim=0)
            return (x * att_weights).sum(dim=0, keepdim=True)
        else:
            raise ValueError(f"Unknown pooling method: {self.pooling}")
    
    def _pool_batch(self, x, batch):
        """배치에 대한 집계"""
        from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool
        
        if self.pooling == 'mean':
            return global_mean_pool(x, batch)
        elif self.pooling == 'sum':
            return global_add_pool(x, batch)
        elif self.pooling == 'max':
            return global_max_pool(x, batch)
        elif self.pooling == 'attention':
            # 배치 처리용 attention pooling
            graph_embeddings = []
            for graph_id in torch.unique(batch):
                mask = (batch == graph_id)
                graph_x = x[mask]
                att_weights = self.attention(graph_x)
                att_weights = F.softmax(att_weights, dim=0)
                graph_emb = (graph_x * att_weights).sum(dim=0, keepdim=True)
                graph_embeddings.append(graph_emb)
            return torch.cat(graph_embeddings, dim=0)
        else:
            raise ValueError(f"Unknown pooling method: {self.pooling}")


def normalize_adjacency(adj):
    """
    인접 행렬을 정규화합니다 (D^(-1/2) A D^(-1/2))
    
    Args:
        adj: 인접 행렬 [N, N]
    
    Returns:
        정규화된 인접 행렬 [N, N]
    """
    # 차수 행렬 계산
    degree = torch.sum(adj, dim=1)
    # 차수의 역수 제곱근 (0으로 나누기 방지)
    degree_inv_sqrt = torch.pow(degree, -0.5)
    degree_inv_sqrt[degree_inv_sqrt == float('inf')] = 0.0
    
    # 대각 행렬 생성
    degree_matrix = torch.diag(degree_inv_sqrt)
    
    # D^(-1/2) A D^(-1/2) 계산
    normalized_adj = torch.mm(torch.mm(degree_matrix, adj), degree_matrix)
    
    return normalized_adj


def create_simple_graph():
    """
    간단한 예제 그래프 생성
    
    Returns:
        Data 객체 (PyTorch Geometric 형식)
    """
    # 노드 특징 (5개 노드, 각각 3차원 특징)
    x = torch.tensor([
        [1.0, 0.0, 0.0],  # 노드 0
        [0.0, 1.0, 0.0],  # 노드 1
        [0.0, 0.0, 1.0],  # 노드 2
        [1.0, 1.0, 0.0],  # 노드 3
        [0.0, 1.0, 1.0],  # 노드 4
    ], dtype=torch.float)
    
    # 엣지 연결 (COO 형식: [2, num_edges])
    # 0-1, 1-2, 2-3, 3-4, 0-3 연결
    edge_index = torch.tensor([
        [0, 1, 1, 2, 2, 3, 3, 4, 0, 3],
        [1, 0, 2, 1, 3, 2, 4, 3, 3, 0]
    ], dtype=torch.long)
    
    # 노드 레이블 (2개 클래스)
    y = torch.tensor([0, 0, 1, 1, 1], dtype=torch.long)
    
    return Data(x=x, edge_index=edge_index, y=y)


def edge_index_to_adjacency(edge_index, num_nodes):
    """
    엣지 인덱스를 인접 행렬로 변환
    
    Args:
        edge_index: 엣지 인덱스 [2, num_edges]
        num_nodes: 노드 개수
    
    Returns:
        인접 행렬 [num_nodes, num_nodes]
    """
    adj = torch.zeros(num_nodes, num_nodes)
    adj[edge_index[0], edge_index[1]] = 1.0
    # 자기 루프 추가
    adj += torch.eye(num_nodes)
    return adj


def get_atomic_features(atomic_number: int, feature_type: str = 'onehot') -> np.ndarray:
    """
    원자 번호로부터 특징 벡터 생성
    
    Args:
        atomic_number: 원자 번호 (1-118)
        feature_type: 특징 타입 ('onehot', 'embedding', 'properties')
    
    Returns:
        특징 벡터
    """
    if feature_type == 'onehot':
        # 간단한 원자 번호 원-핫 인코딩 (주요 원소만)
        common_elements = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                          19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
        if atomic_number in common_elements:
            features = np.zeros(len(common_elements))
            features[common_elements.index(atomic_number)] = 1.0
        else:
            # 알려지지 않은 원소는 마지막 인덱스 사용
            features = np.zeros(len(common_elements))
            features[-1] = 1.0
        return features
    
    elif feature_type == 'embedding':
        # 원자 번호를 정규화된 값으로 사용
        return np.array([atomic_number / 118.0])
    
    elif feature_type == 'properties':
        # 원소의 기본 특성 (주기율표 위치 등)
        # 간단한 버전: 원자 번호, 주기, 족
        period = get_period(atomic_number)
        group = get_group(atomic_number)
        return np.array([atomic_number / 118.0, period / 7.0, group / 18.0])
    
    else:
        raise ValueError(f"Unknown feature_type: {feature_type}")


def get_period(atomic_number: int) -> int:
    """원자 번호로부터 주기 반환"""
    if atomic_number <= 2:
        return 1
    elif atomic_number <= 10:
        return 2
    elif atomic_number <= 18:
        return 3
    elif atomic_number <= 36:
        return 4
    elif atomic_number <= 54:
        return 5
    elif atomic_number <= 86:
        return 6
    else:
        return 7


def get_group(atomic_number: int) -> int:
    """원자 번호로부터 족 반환 (간단한 근사)"""
    # 주기율표의 족 정보 (간단한 근사)
    group_map = {
        1: 1, 2: 2, 3: 1, 4: 2, 5: 13, 6: 14, 7: 15, 8: 16, 9: 17, 10: 18,
        11: 1, 12: 2, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18
    }
    if atomic_number in group_map:
        return group_map[atomic_number]
    # 간단한 근사
    return (atomic_number % 18) + 1


def cif_to_graph(cif_path: str, 
                 cutoff: float = 5.0,
                 feature_type: str = 'onehot',
                 use_periodic: bool = True,
                 add_edge_attr: bool = True,
                 extra_features: Optional[dict] = None) -> Data:
    """
    CIF 파일을 읽어서 그래프로 변환
    
    Args:
        cif_path: CIF 파일 경로
        cutoff: 엣지를 생성할 최대 거리 (Angstrom)
        feature_type: 노드 특징 타입 ('onehot', 'embedding', 'properties')
        use_periodic: 주기적 경계 조건 사용 여부
        add_edge_attr: 엣지 속성(거리) 추가 여부
        extra_features: 추가 노드 특징 딕셔너리 (예: {'bader_charge': [...], 'magmom': [...]})
    
    Returns:
        PyTorch Geometric Data 객체
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE가 설치되지 않았습니다. 설치: pip install ase")
    
    # CIF 파일 읽기
    atoms = ase_read(cif_path)
    
    # ASE arrays에서 추가 특징 추출
    if extra_features is None:
        extra_features = {}
    
    # atoms.arrays에서 추가 정보 추출 (Bader charge, magnetic moment 등)
    for key in atoms.arrays.keys():
        if key not in ['numbers', 'positions']:  # 기본 정보 제외
            extra_features[key] = atoms.arrays[key]
    
    if use_periodic and atoms.pbc.any():
        # 주기적 경계 조건 사용
        return _atoms_to_graph_periodic(atoms, cutoff, feature_type, add_edge_attr, extra_features)
    else:
        # 비주기적 (분자 구조 등)
        return _atoms_to_graph_nonperiodic(atoms, cutoff, feature_type, add_edge_attr, extra_features)


def _atoms_to_graph_periodic(atoms: Atoms, 
                              cutoff: float,
                              feature_type: str,
                              add_edge_attr: bool,
                              extra_features: Optional[dict] = None) -> Data:
    """주기적 경계 조건을 고려한 그래프 생성"""
    positions = atoms.positions
    atomic_numbers = atoms.get_atomic_numbers()
    cell = atoms.cell
    
    num_atoms = len(atoms)
    
    # 노드 특징 생성
    node_features = []
    for i, z in enumerate(atomic_numbers):
        feat = get_atomic_features(z, feature_type)
        
        # 추가 특징 추가 (Bader charge 등)
        if extra_features is not None:
            for key, values in extra_features.items():
                if isinstance(values, (list, np.ndarray)) and len(values) == num_atoms:
                    feat = np.concatenate([feat, [values[i]]])
                elif isinstance(values, dict) and i in values:
                    feat = np.concatenate([feat, [values[i]]])
        
        node_features.append(feat)
    
    x = torch.tensor(np.array(node_features), dtype=torch.float)
    
    # 엣지 생성 (주기적 경계 조건 고려)
    edge_list = []
    edge_attr_list = []
    
    # 모든 원자 쌍에 대해 거리 계산
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            # 주기적 경계 조건을 고려한 최단 거리 계산
            pos_i = positions[i]
            pos_j = positions[j]
            
            # 최소 이미지 컨벤션 (Minimum Image Convention)
            if cell.any():
                # 상대 위치
                rel_pos = pos_j - pos_i
                # 셀 벡터로 변환
                rel_pos_frac = np.linalg.solve(cell.T, rel_pos)
                # [-0.5, 0.5] 범위로 제한
                rel_pos_frac = rel_pos_frac - np.round(rel_pos_frac)
                # 다시 직교 좌표로 변환
                rel_pos = cell.T @ rel_pos_frac
            else:
                rel_pos = pos_j - pos_i
            
            distance = np.linalg.norm(rel_pos)
            
            if distance <= cutoff:
                # 양방향 엣지 추가
                edge_list.append([i, j])
                edge_list.append([j, i])
                
                if add_edge_attr:
                    # 거리를 엣지 속성으로 추가
                    edge_attr_list.append([distance])
                    edge_attr_list.append([distance])
    
    if len(edge_list) == 0:
        # 엣지가 없는 경우 자기 루프만 추가
        edge_list = [[i, i] for i in range(num_atoms)]
        if add_edge_attr:
            edge_attr_list = [[0.0] for _ in range(num_atoms)]
    
    edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
    
    data = Data(x=x, edge_index=edge_index)
    
    if add_edge_attr and len(edge_attr_list) > 0:
        data.edge_attr = torch.tensor(edge_attr_list, dtype=torch.float)
    
    # 추가 정보 저장
    data.positions = torch.tensor(positions, dtype=torch.float)
    data.atomic_numbers = torch.tensor(atomic_numbers, dtype=torch.long)
    
    return data


def _atoms_to_graph_nonperiodic(atoms: Atoms,
                                cutoff: float,
                                feature_type: str,
                                add_edge_attr: bool,
                                extra_features: Optional[dict] = None) -> Data:
    """비주기적 구조를 그래프로 변환"""
    positions = atoms.positions
    atomic_numbers = atoms.get_atomic_numbers()
    
    num_atoms = len(atoms)
    
    # 노드 특징 생성
    node_features = []
    for i, z in enumerate(atomic_numbers):
        feat = get_atomic_features(z, feature_type)
        
        # 추가 특징 추가 (Bader charge 등)
        if extra_features is not None:
            for key, values in extra_features.items():
                if isinstance(values, (list, np.ndarray)) and len(values) == num_atoms:
                    feat = np.concatenate([feat, [values[i]]])
                elif isinstance(values, dict) and i in values:
                    feat = np.concatenate([feat, [values[i]]])
        
        node_features.append(feat)
    
    x = torch.tensor(np.array(node_features), dtype=torch.float)
    
    # 엣지 생성
    edge_list = []
    edge_attr_list = []
    
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = np.linalg.norm(positions[j] - positions[i])
            
            if distance <= cutoff:
                # 양방향 엣지 추가
                edge_list.append([i, j])
                edge_list.append([j, i])
                
                if add_edge_attr:
                    edge_attr_list.append([distance])
                    edge_attr_list.append([distance])
    
    if len(edge_list) == 0:
        # 엣지가 없는 경우 자기 루프만 추가
        edge_list = [[i, i] for i in range(num_atoms)]
        if add_edge_attr:
            edge_attr_list = [[0.0] for _ in range(num_atoms)]
    
    edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
    
    data = Data(x=x, edge_index=edge_index)
    
    if add_edge_attr and len(edge_attr_list) > 0:
        data.edge_attr = torch.tensor(edge_attr_list, dtype=torch.float)
    
    # 추가 정보 저장
    data.positions = torch.tensor(positions, dtype=torch.float)
    data.atomic_numbers = torch.tensor(atomic_numbers, dtype=torch.long)
    
    return data


def example_cif_usage():
    """
    CIF 파일을 사용한 예제
    """
    print("\n" + "=" * 60)
    print("CIF 파일을 그래프로 변환하는 예제")
    print("=" * 60)
    
    if not ASE_AVAILABLE:
        print("\nASE가 설치되지 않았습니다.")
        print("설치: pip install ase")
        return
    
    # 간단한 예제: CIF 파일이 있는 경우
    print("\n사용법:")
    print("  from gnn import cif_to_graph")
    print("  data = cif_to_graph('structure.cif', cutoff=5.0)")
    print("  print(f'노드 개수: {data.x.shape[0]}')")
    print("  print(f'엣지 개수: {data.edge_index.shape[1] // 2}')")
    
    print("\n주요 파라미터:")
    print("  - cutoff: 엣지를 생성할 최대 거리 (기본값: 5.0 Angstrom)")
    print("  - feature_type: 노드 특징 타입 ('onehot', 'embedding', 'properties')")
    print("  - use_periodic: 주기적 경계 조건 사용 여부 (기본값: True)")
    print("  - add_edge_attr: 엣지 속성(거리) 추가 여부 (기본값: True)")
    
    return


def example_json_usage():
    """
    JSON 파일을 사용한 예제
    """
    print("\n" + "=" * 60)
    print("JSON 파일을 그래프로 변환하는 예제")
    print("=" * 60)
    
    if not ASE_AVAILABLE:
        print("\nASE가 설치되지 않았습니다.")
        print("설치: pip install ase")
        return
    
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
    print("  print(f'노드 특징 차원: {data.x.shape[1]}')  # 원소 특징 + 추가 특징")
    print("  print(f'엣지 개수: {data.edge_index.shape[1] // 2}')")
    
    print("\nJSON 파일 형식 예시:")
    print("  {")
    print("    \"atoms\": {...},  # ASE 형식 구조")
    print("    \"extra\": {")
    print("      \"bader_charge\": [0.1, 0.2, -0.3, ...],  # 각 원자별 Bader charge")
    print("      \"magmom\": [1.0, -1.0, 0.0, ...]  # 각 원자별 자기 모멘트")
    print("    }")
    print("  }")
    
    print("\n또는 ASE arrays에 저장된 경우:")
    print("  atoms.arrays['bader_charge'] = [0.1, 0.2, ...]")
    print("  atoms.write('structure.json')")
    
    print("\n주요 파라미터:")
    print("  - cutoff: 엣지를 생성할 최대 거리 (기본값: 5.0 Angstrom)")
    print("  - feature_type: 노드 특징 타입 ('onehot', 'embedding', 'properties')")
    print("  - use_periodic: 주기적 경계 조건 사용 여부 (기본값: True)")
    print("  - extra_feature_keys: 추출할 추가 특징 키 리스트")
    print("    (None이면 자동으로 모든 가능한 특징 추출)")
    
    return


def example_cgnn_usage():
    """
    Crystal Graph Neural Network 사용 예제
    """
    print("\n" + "=" * 60)
    print("Crystal Graph Neural Network (CGNN) 사용 예제")
    print("=" * 60)
    
    print("\nCGNN은 결정 구조의 특성을 예측하기 위한 특화된 모델입니다.")
    print("(예: 형성 에너지, 밴드갭, 탄성 계수 등)")
    
    print("\n사용법:")
    print("  from gnn import CrystalGraphNeuralNetwork, cif_to_graph")
    print("  ")
    print("  # 1. 결정 구조를 그래프로 변환")
    print("  data = cif_to_graph('structure.cif', cutoff=5.0, add_edge_attr=True)")
    print("  ")
    print("  # 2. CGNN 모델 생성")
    print("  model = CrystalGraphNeuralNetwork(")
    print("      num_features=data.x.shape[1],  # 노드 특징 차원")
    print("      hidden_dim=64,")
    print("      num_outputs=1,  # 회귀 작업 (형성 에너지 등)")
    print("      num_layers=3,")
    print("      pooling='mean'  # 또는 'sum', 'max', 'attention'")
    print("  )")
    print("  ")
    print("  # 3. 예측")
    print("  output = model(data.x, data.edge_index, data.edge_attr)")
    print("  print(f'예측값: {output.item()}')")
    
    print("\nCGNN의 주요 특징:")
    print("  - Crystal Graph Convolution: 거리 기반 가중 메시지 전달")
    print("  - 엣지 속성 활용: 원자 간 거리를 엣지 특징으로 사용")
    print("  - 그래프 레벨 집계: 전체 결정 구조의 특성 예측")
    print("  - 주기적 경계 조건 고려: 결정 구조의 특성 반영")
    
    print("\n주요 파라미터:")
    print("  - num_features: 입력 노드 특징 차원")
    print("  - hidden_dim: 은닉층 차원 (일반적으로 64-128)")
    print("  - num_outputs: 출력 차원 (회귀: 1, 분류: 클래스 수)")
    print("  - num_layers: 레이어 개수 (일반적으로 3-5)")
    print("  - pooling: 그래프 집계 방법 ('mean', 'sum', 'max', 'attention')")
    print("  - use_edge_attr: 엣지 속성 사용 여부 (기본값: True)")
    
    return


def json_to_graph(json_path: str,
                  cutoff: float = 5.0,
                  feature_type: str = 'onehot',
                  use_periodic: bool = True,
                  add_edge_attr: bool = True,
                  extra_feature_keys: Optional[List[str]] = None) -> Data:
    """
    JSON 파일을 읽어서 그래프로 변환
    
    JSON 파일 형식:
    - ASE 형식: ASE가 읽을 수 있는 표준 형식
    - 커스텀 형식: 'atoms' 키에 ASE 형식 구조, 'extra' 키에 추가 정보
    
    Args:
        json_path: JSON 파일 경로
        cutoff: 엣지를 생성할 최대 거리 (Angstrom)
        feature_type: 노드 특징 타입 ('onehot', 'embedding', 'properties')
        use_periodic: 주기적 경계 조건 사용 여부
        add_edge_attr: 엣지 속성(거리) 추가 여부
        extra_feature_keys: JSON에서 추출할 추가 특징 키 리스트 
                          (예: ['bader_charge', 'magmom', 'forces'])
                          None이면 자동으로 모든 가능한 특징 추출
    
    Returns:
        PyTorch Geometric Data 객체
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE가 설치되지 않았습니다. 설치: pip install ase")
    
    # JSON 파일 읽기
    with open(json_path, 'r') as f:
        json_data = json.load(f)
    
    # ASE atoms 객체 생성
    # ASE가 직접 JSON을 읽을 수 있음
    atoms = ase_read(json_path)
    
    # 추가 특징 추출
    extra_features = {}
    
    # 1. ASE arrays에서 추출 (Bader charge, magnetic moment 등)
    for key in atoms.arrays.keys():
        if key not in ['numbers', 'positions']:
            if extra_feature_keys is None or key in extra_feature_keys:
                extra_features[key] = atoms.arrays[key]
    
    # 2. JSON 파일에서 직접 추출 (커스텀 형식)
    if isinstance(json_data, dict):
        # 'extra' 키가 있는 경우
        if 'extra' in json_data and isinstance(json_data['extra'], dict):
            for key, value in json_data['extra'].items():
                if extra_feature_keys is None or key in extra_feature_keys:
                    if isinstance(value, list) and len(value) == len(atoms):
                        extra_features[key] = np.array(value)
        
        # 최상위 레벨에 원자별 정보가 있는 경우
        # 예: {'bader_charge': [0.1, 0.2, ...], 'magmom': [1.0, -1.0, ...]}
        for key, value in json_data.items():
            if key not in ['1', 'ids', 'nextid']:  # ASE 메타데이터 제외
                if isinstance(value, list) and len(value) == len(atoms):
                    if extra_feature_keys is None or key in extra_feature_keys:
                        extra_features[key] = np.array(value)
                elif isinstance(value, dict):
                    # 인덱스 기반 딕셔너리인 경우
                    if all(isinstance(k, (int, str)) and k.isdigit() if isinstance(k, str) else True 
                           for k in value.keys()):
                        if extra_feature_keys is None or key in extra_feature_keys:
                            extra_features[key] = value
    
    # 그래프로 변환
    if use_periodic and atoms.pbc.any():
        return _atoms_to_graph_periodic(atoms, cutoff, feature_type, add_edge_attr, extra_features)
    else:
        return _atoms_to_graph_nonperiodic(atoms, cutoff, feature_type, add_edge_attr, extra_features)


def batch_json_to_graphs(json_paths: List[str],
                         cutoff: float = 5.0,
                         feature_type: str = 'onehot',
                         use_periodic: bool = True,
                         extra_feature_keys: Optional[List[str]] = None) -> List[Data]:
    """
    여러 JSON 파일을 일괄 처리하여 그래프 리스트로 변환
    
    Args:
        json_paths: JSON 파일 경로 리스트
        cutoff: 엣지를 생성할 최대 거리
        feature_type: 노드 특징 타입
        use_periodic: 주기적 경계 조건 사용 여부
        extra_feature_keys: 추출할 추가 특징 키 리스트
    
    Returns:
        Data 객체 리스트
    """
    graphs = []
    for json_path in json_paths:
        try:
            graph = json_to_graph(json_path, cutoff, feature_type, use_periodic, 
                                 extra_feature_keys=extra_feature_keys)
            graphs.append(graph)
        except Exception as e:
            print(f"Warning: {json_path} 처리 중 오류 발생: {e}")
            continue
    return graphs


def batch_cif_to_graphs(cif_paths: List[str], 
                        cutoff: float = 5.0,
                        feature_type: str = 'onehot',
                        use_periodic: bool = True) -> List[Data]:
    """
    여러 CIF 파일을 일괄 처리하여 그래프 리스트로 변환
    
    Args:
        cif_paths: CIF 파일 경로 리스트
        cutoff: 엣지를 생성할 최대 거리
        feature_type: 노드 특징 타입
        use_periodic: 주기적 경계 조건 사용 여부
    
    Returns:
        Data 객체 리스트
    """
    graphs = []
    for cif_path in cif_paths:
        try:
            graph = cif_to_graph(cif_path, cutoff, feature_type, use_periodic)
            graphs.append(graph)
        except Exception as e:
            print(f"Warning: {cif_path} 처리 중 오류 발생: {e}")
            continue
    return graphs


def example_basic_usage():
    """
    기본 사용 예제
    """
    print("=" * 60)
    print("GNN 기본 사용 예제")
    print("=" * 60)
    
    # 그래프 생성
    data = create_simple_graph()
    num_nodes = data.x.shape[0]
    num_features = data.x.shape[1]
    num_classes = len(torch.unique(data.y))
    
    print(f"\n그래프 정보:")
    print(f"  - 노드 개수: {num_nodes}")
    print(f"  - 특징 차원: {num_features}")
    print(f"  - 클래스 개수: {num_classes}")
    print(f"  - 엣지 개수: {data.edge_index.shape[1] // 2}")
    
    # 인접 행렬 생성 및 정규화
    adj = edge_index_to_adjacency(data.edge_index, num_nodes)
    normalized_adj = normalize_adjacency(adj)
    
    # 모델 생성
    hidden_dim = 16
    model = SimpleGCN(num_features, hidden_dim, num_classes)
    
    print(f"\n모델 구조:")
    print(f"  - 입력 차원: {num_features}")
    print(f"  - 은닉 차원: {hidden_dim}")
    print(f"  - 출력 차원: {num_classes}")
    
    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(data.x, normalized_adj)
        predictions = output.argmax(dim=1)
    
    print(f"\n예측 결과:")
    print(f"  실제 레이블: {data.y.tolist()}")
    print(f"  예측 레이블: {predictions.tolist()}")
    print(f"  예측 확률:\n{torch.exp(output).numpy()}")
    
    return model, data, normalized_adj


def example_with_pyg():
    """
    PyTorch Geometric을 사용한 고급 예제
    """
    print("\n" + "=" * 60)
    print("PyTorch Geometric 사용 예제")
    print("=" * 60)
    
    # 간단한 그래프 생성
    data = create_simple_graph()
    
    # PyTorch Geometric의 GCN 레이어 사용
    class PyGGCN(nn.Module):
        def __init__(self, num_features, hidden_dim, num_classes):
            super(PyGGCN, self).__init__()
            self.conv1 = GCNConv(num_features, hidden_dim)
            self.conv2 = GCNConv(hidden_dim, num_classes)
        
        def forward(self, x, edge_index):
            x = self.conv1(x, edge_index)
            x = F.relu(x)
            x = F.dropout(x, training=self.training)
            x = self.conv2(x, edge_index)
            return F.log_softmax(x, dim=1)
    
    model = PyGGCN(data.x.shape[1], 16, len(torch.unique(data.y)))
    
    print(f"\nPyG GCN 모델 생성 완료")
    print(f"  - 파라미터 개수: {sum(p.numel() for p in model.parameters())}")
    
    # Forward pass
    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index)
        predictions = output.argmax(dim=1)
    
    print(f"\n예측 결과:")
    print(f"  실제 레이블: {data.y.tolist()}")
    print(f"  예측 레이블: {predictions.tolist()}")
    
    return model


def example_training():
    """
    간단한 학습 예제
    """
    print("\n" + "=" * 60)
    print("GNN 학습 예제")
    print("=" * 60)
    
    # 데이터 준비
    data = create_simple_graph()
    num_nodes = data.x.shape[0]
    adj = edge_index_to_adjacency(data.edge_index, num_nodes)
    normalized_adj = normalize_adjacency(adj)
    
    # 모델 및 옵티마이저
    model = SimpleGCN(data.x.shape[1], 16, len(torch.unique(data.y)))
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = nn.NLLLoss()
    
    # 학습
    model.train()
    num_epochs = 100
    
    print(f"\n학습 시작 (epochs: {num_epochs})...")
    for epoch in range(num_epochs):
        optimizer.zero_grad()
        output = model(data.x, normalized_adj)
        loss = criterion(output, data.y)
        loss.backward()
        optimizer.step()
        
        if (epoch + 1) % 20 == 0:
            accuracy = (output.argmax(dim=1) == data.y).float().mean()
            print(f"  Epoch {epoch+1:3d}: Loss = {loss.item():.4f}, Accuracy = {accuracy:.4f}")
    
    # 평가
    model.eval()
    with torch.no_grad():
        output = model(data.x, normalized_adj)
        predictions = output.argmax(dim=1)
        accuracy = (predictions == data.y).float().mean()
    
    print(f"\n최종 결과:")
    print(f"  정확도: {accuracy:.4f}")
    print(f"  실제 레이블: {data.y.tolist()}")
    print(f"  예측 레이블: {predictions.tolist()}")


if __name__ == "__main__":
    # 기본 사용 예제
    example_basic_usage()
    
    # PyTorch Geometric 예제
    try:
        example_with_pyg()
    except ImportError:
        print("\nPyTorch Geometric이 설치되지 않았습니다.")
        print("설치: pip install torch-geometric")
    
    # 학습 예제
    example_training()
    
    # CIF 파일 사용 예제
    example_cif_usage()
    
    # JSON 파일 사용 예제
    example_json_usage()
    
    # CGNN 사용 예제
    example_cgnn_usage()
    
    print("\n" + "=" * 60)
    print("예제 실행 완료!")
    print("=" * 60)
