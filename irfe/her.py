import pandas as pd
import matplotlib.pyplot as plt
import os

# 기본 경로 설정 (hbe.py, oer.py와 동일)
base_path = '/Users/hailey/Desktop/4_IrFe3'
figures_dir = os.path.join(base_path, 'figures')
os.makedirs(figures_dir, exist_ok=True)

# H binding energy 데이터 로드
hbe_csv = os.path.join(base_path, 'H_binding_energies_final.csv')
df = pd.read_csv(hbe_csv)

# 필요한 값 추출: IrFe (low), Ir (low), Ir (high)
# low: 4_atom_top2 (IrFe), 5_atom_top (Ir), high: 1_layer_top (Ir)
irfe_low = df[(df['System_name'] == 'IrFe') & (df['Site'] == '4_atom_top2')]['Binding_energy'].values[0]
ir_low = df[(df['System_name'] == 'Ir') & (df['Site'] == '5_atom_top')]['Binding_energy'].values[0]
ir_high = df[(df['System_name'] == 'Ir') & (df['Site'] == '1_layer_top')]['Binding_energy'].values[0]

# HER energetics: 2H+ + 2e- -> H2 (ΔG_H for each step)
# 일반적으로 HER는 ΔG_H (H* -> 0.5H2)로 표현
# 여기서는 각 표면의 ΔG_H를 1-step reaction profile로 그림
her_labels = ['Ir (low coverage, η = 0.26 V)', 'Ir (high coverage, η = 0.21 V)', 'IrFe (low coverage, η = 0.06 V)']
her_dg = [ir_low, ir_high, irfe_low]
colors = ['black', 'black', 'red']
linestyles = ['--', '-', '-']
zorders = [2, 2, 1]

# 꺾은선 그래프 (step diagram) 그리기
def plot_her_profile(her_dg, her_labels, figures_dir):
    plt.figure(figsize=(4, 3))
    for i, (dg, label) in enumerate(zip(her_dg, her_labels)):
        # HER: 0 -> ΔG_H -> 0 (H2)
        yy = [0, 0, dg, dg, 0, 0]
        xx = [0, 1, 2, 3, 4, 5]
        plt.plot(xx, yy, marker='', label=label, zorder=zorders[i], color=colors[i], linestyle=linestyles[i])
    # plt.axhline(y=0, color='black', linestyle='--', linewidth=1.0, zorder=1)
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('Relative Energy (ΔG, eV)')
    plt.xticks([0.5, 2.5, 4.5], ['*', '*H', '1/2 H$_2$'])
    plt.xlim(-0.2, 5.2)
    plt.ylim(-0.6, 0.1)
    plt.legend(loc='lower center', frameon=True, framealpha=1.0, ncol=1, handletextpad=0.5, labelspacing=0.2)
    plt.tight_layout()
    plot_path = os.path.join(figures_dir, 'HER_profile.png')
    plt.savefig(plot_path, bbox_inches='tight', transparent=True, dpi=300)
    print(f"HER profile plot saved to {plot_path}")
    plt.close()

if __name__ == "__main__":
    plot_her_profile(her_dg, her_labels, figures_dir) 