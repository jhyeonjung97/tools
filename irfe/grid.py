from ase.io import read
import numpy as np

atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()

# 음의 방향으로 셀 크기 변경
for i in range(5):
    scale_x = 1.0 - (i + 1) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x
    
    for j in range(5):
        scale_y = 1.0 - (j + 1) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y
        
        # 원자 위치를 비율에 맞게 스케일링
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y
        
        atoms.write(f'n{i}{j}.vasp')

# 원본 상태로 복원
atoms = read('CONTCAR')
original_cell = atoms.cell.copy()
original_positions = atoms.positions.copy()

# 양의 방향으로 셀 크기 변경
for i in range(5):
    scale_x = 1.0 + (i + 1) * 0.1 / np.linalg.norm(original_cell[0])
    atoms.cell[0] = original_cell[0] * scale_x
    
    for j in range(5):
        scale_y = 1.0 + (j + 1) * 0.1 / np.linalg.norm(original_cell[1])
        atoms.cell[1] = original_cell[1] * scale_y
        
        # 원자 위치를 비율에 맞게 스케일링
        atoms.positions[:, 0] = original_positions[:, 0] * scale_x
        atoms.positions[:, 1] = original_positions[:, 1] * scale_y
        
        atoms.write(f'p{i}{j}.vasp')

# mkdir n00 n01 n02 n03 n04
# mkdir n10 n11 n12 n13 n14
# mkdir n20 n21 n22 n23 n24
# mkdir n30 n31 n32 n33 n34
# mkdir n40 n41 n42 n43 n44

# mkdir p00 p01 p02 p03 p04
# mkdir p10 p11 p12 p13 p14
# mkdir p20 p21 p22 p23 p24
# mkdir p30 p31 p32 p33 p34
# mkdir p40 p41 p42 p43 p44

