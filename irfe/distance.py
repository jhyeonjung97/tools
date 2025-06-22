from ase.io import read
import numpy as np

atoms = read('CONTCAR')

# Ir 원자들의 인덱스 찾기
ir_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'Ir']

if ir_indices:
    # Ir 원자들의 z 좌표 가져오기
    ir_positions = atoms.get_positions()[ir_indices]
    z_coords = ir_positions[:, 2]  # z 좌표 (3번째 열)
    
    # 가장 높은 z 좌표를 가진 Ir 원자 찾기
    max_z_index = np.argmax(z_coords)
    highest_ir_index = ir_indices[max_z_index]
    highest_z = z_coords[max_z_index]
    
    # print(f"가장 높은 Ir 원자:")
    # print(f"  인덱스: {highest_ir_index}")
    # print(f"  위치: {atoms[highest_ir_index].position}")
    # print(f"  z 좌표: {highest_z:.6f} Å")
    
    # 금속 원자들(Ir, Fe)의 인덱스 찾기
    metal_indices = [i for i, atom in enumerate(atoms) if atom.symbol in ['Ir', 'Fe']]
    
    # 기준 Ir 원자의 위치
    reference_position = atoms[highest_ir_index].position
    
    # 모든 금속 원자와의 거리 계산
    distances = []
    for idx in metal_indices:
        if idx != highest_ir_index:  # 자기 자신 제외
            distance = np.linalg.norm(atoms[idx].position - reference_position)
            distances.append((idx, distance, atoms[idx].symbol))
    
    # 거리순으로 정렬
    distances.sort(key=lambda x: x[1])
    
    # 가장 가까운 3개 선택
    closest_3 = distances[:3]
    
    # print(f"\n가장 가까운 금속 원자 3개:")
    # for i, (idx, distance, symbol) in enumerate(closest_3, 1):
    #     print(f"  {i}번째: {symbol} 원자 (인덱스: {idx})")
    #     print(f"    위치: {atoms[idx].position}")
    #     print(f"    거리: {distance:.6f} Å")
    
    # 거리들의 평균 계산
    avg_distance = np.mean([dist[1] for dist in closest_3])
    # print(f"\n가장 가까운 금속 원자 3개와의 평균 거리: {avg_distance:.6f} Å")
    print(avg_distance)
else:
    print("CONTCAR 파일에서 Ir 원자를 찾을 수 없습니다.")






