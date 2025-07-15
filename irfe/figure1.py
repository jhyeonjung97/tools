import os
import pandas as pd
import matplotlib.pyplot as plt

# 파일 경로 설정
base_path = "/Users/hailey/Desktop/4_IrFe3"
ir_file = os.path.join(base_path, "pourbaix_data_Ir.csv")
irfe_file = os.path.join(base_path, "pourbaix_data_IrFe.csv")

# 데이터 읽기
ir = pd.read_csv(ir_file)
irfe = pd.read_csv(irfe_file)

# 색상 설정
color_dict = {
    'H': 'blue',
    'OH': 'green',
    'O': 'red'
}

# 전압 범위
Umin, Umax = -0.5, 2.0
U_range = [Umin, Umax]

plt.figure(figsize=(4.5, 3))

# IrFe: 실선
for ads in ['H', 'OH', 'O']:
    row = irfe[irfe['ads'] == ads].iloc[0]
    plt.plot(U_range, [row['Gmin'], row['Gmax']], color=color_dict[ads], 
    label=f"IrFe*{ads}", linestyle='-')

# x=0 수직선
plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
# clean 기준선
plt.plot(U_range, [0, 0], color='black', label='IrFe, Ir clean')

# Ir: 점선
for ads in ['H', 'OH', 'O']:
    row = ir[ir['ads'] == ads].iloc[0]
    plt.plot(U_range, [row['Gmin'], row['Gmax']], color=color_dict[ads], 
    label=f"Ir*{ads}", linestyle='--', linewidth=1.0)

plt.xlabel('Potential (V vs RHE)')
plt.ylabel('ΔG (eV/site)')
plt.xlim(Umin, Umax)
plt.ylim(-1.5, 0.5)
leg = plt.legend(
    loc='lower left',
    bbox_to_anchor=(0.0, 0.0),
    handletextpad=0.3,
    labelspacing=0.2,
    ncol=2,
    columnspacing=-1.8,
    frameon=True,
    framealpha=1.0,
    handlelength=1.4  # 기본값은 2.0, 더 작게 하면 선이 짧아집니다
)
plt.tight_layout()

# 파일로 저장
output_path = os.path.join(base_path, "figure1.png")
plt.savefig(output_path, bbox_inches='tight', dpi=300, transparent=True)
plt.close()
print(f"Plot saved to {output_path}") 