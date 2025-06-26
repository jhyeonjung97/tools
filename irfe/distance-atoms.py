from ase.io import read
import numpy as np
from itertools import combinations

atoms = read('CONTCAR')

# 총 원자 개수 출력
total_atoms = len(atoms)
print(f"총 원자 개수: {total_atoms}")

# periodic 조건 확인
pbc = atoms.get_pbc()
print(f"Periodic 조건: {pbc}")

# z축 방향으로 낮은 24개 원자 선정
z_positions = atoms.get_positions()[:, 2]  # z 좌표만 추출
lowest_indices = np.argsort(z_positions)[:24]  # z값이 가장 낮은 24개 원자의 인덱스

print(f"선정된 원자 개수: {len(lowest_indices)}")

# 선정된 원자들 간의 거리 계산 (periodic 조건 고려)
distances = []
for i, j in combinations(lowest_indices, 2):
    # periodic=True로 명시적으로 설정하여 periodic 조건 고려
    distance = atoms.get_distance(i, j, mic=True)  # mic=True는 "minimum image convention"을 의미
    distances.append(distance)

# 3 옴스트롱 이하의 거리만 필터링
distances_filtered = [d for d in distances if d <= 3.0]

# 평균 거리 계산
average_distance = np.mean(distances_filtered)

print(f"z축 방향으로 낮은 24개 원자 간의 평균 거리: {average_distance:.4f} Å")
print(f"총 거리 쌍의 수: {len(distances)}")
print(f"3 Å 이하 거리 쌍의 수: {len(distances_filtered)}")
print(f"3 Å 이하 거리들의 평균: {average_distance:.4f} Å")
print(f"최소 거리: {min(distances):.4f} Å")
print(f"최대 거리: {max(distances):.4f} Å")
print(f"3 Å 이하 최소 거리: {min(distances_filtered):.4f} Å")
print(f"3 Å 이하 최대 거리: {max(distances_filtered):.4f} Å")



