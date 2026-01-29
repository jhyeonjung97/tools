# Graph Neural Network (GNN) 구현 예제

이 디렉토리는 Graph Neural Network의 기본 구현 예제를 포함합니다.

## 설치

```bash
pip install -r requirements.txt
```

## 주요 내용

### 1. BasicGCNLayer
- 기본적인 Graph Convolutional Network 레이어 구현
- 수식: H^(l+1) = σ(D^(-1/2) A D^(-1/2) H^(l) W^(l))

### 2. SimpleGCN
- 2-layer GCN 모델
- 노드 분류 작업에 사용 가능

### 3. 유틸리티 함수
- `normalize_adjacency`: 인접 행렬 정규화
- `edge_index_to_adjacency`: 엣지 인덱스를 인접 행렬로 변환
- `create_simple_graph`: 예제 그래프 생성

## 사용 예제

### 기본 실행
```bash
python gnn.py
```

### JSON 파일을 그래프로 변환 (권장)
**각 원자별 정보를 포함한 CGNN 사용:**
```python
from gnn import json_to_graph, CrystalGraphNeuralNetwork

# JSON 파일을 그래프로 변환 (Bader charge 등 포함)
data = json_to_graph('structure.json', cutoff=5.0,
                     extra_feature_keys=['bader_charge', 'magmom', 'coordination'])

print(f"노드 개수: {data.x.shape[0]}")
print(f"노드 특징 차원: {data.x.shape[1]}")  # 원소 특징 + 추가 특징들
print(f"엣지 개수: {data.edge_index.shape[1] // 2}")

# CGNN 모델 생성 및 예측
model = CrystalGraphNeuralNetwork(
    num_features=data.x.shape[1],  # 원소 특징 + 추가 특징들
    hidden_dim=64,
    num_outputs=1  # 형성 에너지 등 예측
)
output = model(data.x, data.edge_index, data.edge_attr)
```

### CIF 파일을 그래프로 변환
```python
from gnn import cif_to_graph

# CIF 파일을 그래프로 변환
data = cif_to_graph('structure.cif', cutoff=5.0)

# 추가 정보를 넣으려면 extra_features 사용
extra_features = {
    'bader_charge': [0.1, 0.2, -0.3, ...],  # 각 원자별 정보
    'magmom': [1.0, -1.0, 0.0, ...]
}
data = cif_to_graph('structure.cif', cutoff=5.0, extra_features=extra_features)
```

### 여러 JSON 파일 일괄 처리
```python
from gnn import batch_json_to_graphs

json_files = ['file1.json', 'file2.json', 'file3.json']
graphs = batch_json_to_graphs(json_files, cutoff=5.0,
                              extra_feature_keys=['bader_charge', 'magmom'])
```

### 코드에서 사용
```python
from gnn import SimpleGCN, create_simple_graph, normalize_adjacency, edge_index_to_adjacency

# 데이터 준비
data = create_simple_graph()
adj = edge_index_to_adjacency(data.edge_index, data.x.shape[0])
normalized_adj = normalize_adjacency(adj)

# 모델 생성 및 사용
model = SimpleGCN(num_features=3, hidden_dim=16, num_classes=2)
output = model(data.x, normalized_adj)
```

### CIF 파일 예제 실행
```bash
python example_cif.py
```

### JSON 파일을 그래프로 변환
```python
from gnn import json_to_graph

# 기본 사용
data = json_to_graph('structure.json', cutoff=5.0)

# Bader charge 등 추가 특징 포함
data = json_to_graph('structure.json', cutoff=5.0,
                     extra_feature_keys=['bader_charge', 'magmom'])

print(f"노드 개수: {data.x.shape[0]}")
print(f"노드 특징 차원: {data.x.shape[1]}")  # 원소 특징 + 추가 특징
print(f"엣지 개수: {data.edge_index.shape[1] // 2}")
```

### 여러 JSON 파일 일괄 처리
```python
from gnn import batch_json_to_graphs

json_files = ['file1.json', 'file2.json', 'file3.json']
graphs = batch_json_to_graphs(json_files, cutoff=5.0,
                              extra_feature_keys=['bader_charge'])
```

### JSON 파일 예제 실행
```bash
python example_json.py
```

### Crystal Graph Neural Network (CGNN) 사용
```python
from gnn import CrystalGraphNeuralNetwork, cif_to_graph

# 결정 구조를 그래프로 변환
data = cif_to_graph('structure.cif', cutoff=5.0, add_edge_attr=True)

# CGNN 모델 생성
model = CrystalGraphNeuralNetwork(
    num_features=data.x.shape[1],
    hidden_dim=64,
    num_outputs=1,  # 회귀 작업 (형성 에너지 등)
    num_layers=3,
    pooling='mean'  # 또는 'sum', 'max', 'attention'
)

# 예측
output = model(data.x, data.edge_index, data.edge_attr)
print(f'예측값: {output.item()}')
```

### CGNN 예제 실행
```bash
python example_cgnn.py
```

### JSON 파일을 사용한 CGNN 예제 (권장)
```bash
python example_cgnn_json.py
```

이 예제는 CGNN의 핵심 개념을 보여줍니다:
- 각 원자별 정보(Bader charge, magnetic moment 등)를 노드 특징으로 사용
- 노드 간 그래프 관계를 이용한 결정 구조 표현

## 주요 기능

### 1. Crystal Graph Neural Network (CGNN)
**CGNN의 핵심 개념:**
- **각 노드(원자)에 대한 정보를 직접 넣을 수 있음**: Bader charge, magnetic moment, coordination number 등 각 원자별 정보를 노드 특징으로 사용
- **노드 간 그래프 관계를 이용**: 원자 간 거리 기반 엣지로 결정 구조의 연결성 표현
- **거리 기반 가중 메시지 전달**: 가까운 원자일수록 더 큰 영향
- **그래프 레벨 집계**: 전체 결정 구조의 특성 예측 (형성 에너지, 밴드갭 등)

**구현 특징:**
- Materials Project 스타일의 CGNN 아키텍처
- Crystal Graph Convolution 레이어
- 다양한 Pooling 방법: mean, sum, max, attention

### 2. JSON 파일 지원 (권장)
**JSON 파일을 사용하는 이유:**
- ✅ **각 원자별 정보 추가가 쉬움**: Bader charge, magnetic moment 등 원자별 데이터를 쉽게 포함 가능
- ✅ **유연한 데이터 구조**: 원하는 모든 정보를 JSON에 저장 가능
- ✅ **ASE 형식 및 커스텀 형식 지원**: 기존 ASE JSON 파일도 그대로 사용 가능

**JSON 파일 형식 예시:**
```json
{
  "atoms": {...},  // ASE 형식 구조
  "extra": {
    "bader_charge": [0.1, 0.2, -0.3, ...],  // 각 원자별 Bader charge
    "magmom": [1.0, -1.0, 0.0, ...],        // 각 원자별 자기 모멘트
    "coordination": [4, 6, 4, ...]          // 각 원자별 배위수
  }
}
```

### 3. CIF 파일 지원
- `.cif` 파일도 지원하지만, 각 원자별 추가 정보를 넣기 위해서는 ASE arrays를 사용해야 함
- JSON 파일이 더 직관적이고 편리함

### 4. 공통 기능
- 주기적 경계 조건 지원
- 원자 간 거리 기반 엣지 생성
- 다양한 노드 특징 타입 지원 (onehot, embedding, properties)

### 2. 그래프 변환
- 원자 = 노드
- 원자 간 결합 = 엣지 (cutoff 거리 기반)
- 노드 특징: 원소 종류, 원소 특성 등
- 엣지 특징: 원자 간 거리 (선택적)

### 3. 예제 그래프 구조

```
0 -- 1 -- 2 -- 3 -- 4
|         |
+---------+
```

- 노드 개수: 5
- 특징 차원: 3
- 클래스 개수: 2

## 참고 자료

- [Graph Convolutional Networks (Kipf & Welling, 2017)](https://arxiv.org/abs/1609.02907)
- [PyTorch Geometric Documentation](https://pytorch-geometric.readthedocs.io/)
