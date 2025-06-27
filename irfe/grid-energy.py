from ase.io import read
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.DataFrame(columns=[-4, -3, -2, -1, 0, 1, 2, 3, 4], index=[4, 3, 2, 1, 0, -1, -2, -3, -4])
for k in range(1, 5):
    for i in range(5):
        for j in range(5):
            try:
                atoms = read(f'{k}_/{i}{j}/final.json')
                energy = atoms.get_potential_energy()
            except:
                energy = None
            if k == 1:
                df.loc[+j, +i] = energy
            elif k == 2:
                df.loc[+j, -i] = energy
            elif k == 3:
                df.loc[-j, -i] = energy
            elif k == 4:
                df.loc[-j, +i] = energy

energy_reference = df.loc[0, 0]
df = df - energy_reference
print(df)
df.to_csv('grid.tsv', sep='\t', float_format='%.2f')

# 히트맵 시각화
plt.figure(figsize=(10, 8))

# 데이터를 numpy 배열로 변환하고 None 값을 NaN으로 처리
data_matrix = df.values.astype(float)  # float 타입으로 변환
data_matrix = np.where(pd.isna(data_matrix), np.nan, data_matrix)  # None을 NaN으로 처리

# 히트맵 생성 (음수는 파랑, 양수는 빨강)
im = plt.imshow(data_matrix, cmap='RdBu_r', aspect='equal', interpolation='nearest', vmin=-1, vmax=1)

# 컬러바 추가
cbar = plt.colorbar(im)
cbar.set_label('Energy (eV)', rotation=270, labelpad=15)

# 축 레이블 설정
plt.xticks(range(len(df.columns)), df.columns)
plt.yticks(range(len(df.index)), df.index)
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.title('Energy Grid Heatmap')

# 각 셀에 값 표시 (NaN이 아닌 경우에만)
for i in range(len(df.columns)):  # y coordinate (row)
    for j in range(len(df.index)):  # x coordinate (column)
        if not np.isnan(data_matrix[j, i]):
            plt.text(i, j, f'{data_matrix[j, i]:.3f}', 
                    ha='center', va='center', fontsize=8, 
                    color='black' if i == 0 and j == 0 else 'white')

plt.tight_layout()
plt.savefig('energy_grid_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()