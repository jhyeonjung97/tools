from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

root = "~/Desktop/3_RuO2/1_RuO2_Ueff"

# gas (전역 상수: 두 스크립트에서 동일 값 사용)
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
go = gh2o - gh2 - 1.229*2

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo


def get_energy_from_json(json_path):
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None


def calculate_o_and_oh_adsorption_energies():
    root_path = Path(root).expanduser()
    v_o_path = root_path / "3_O_V"
    o_ohcus_path = root_path / "4_O_OH"
    o_o_path = root_path / "5_O_O"

    ueff_values = [0, 1, 2, 3, 4]
    energies_v_o = []
    energies_o_ohcus = []
    energies_o_o = []

    # V_O 에너지 추출
    for ueff in ueff_values:
        json_path = v_o_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_v_o.append(energy)
    
    # O_OHcus 에너지 추출
    for ueff in ueff_values:
        json_path = o_ohcus_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_o_ohcus.append(energy)
    
    # O_O 에너지 추출
    for ueff in ueff_values:
        json_path = o_o_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_o_o.append(energy)

    if len(energies_v_o) == len(energies_o_ohcus) == len(energies_o_o) == len(ueff_values):
        # 두 번째 O와 OH 흡착 에너지만 계산
        O_Ocovered = []    # O_O - V_O
        OH_Ocovered = []   # O_OHcus - V_O

        for i, ueff in enumerate(ueff_values):
            # 기준 가스 상
            o_gas_energy = gh2o - gh2
            oh_gas_energy = gh2o - gh2/2

            # 두 번째 O 흡착 에너지
            O_Ocovered.append(energies_o_o[i] + dgo - energies_v_o[i] - o_gas_energy)
            
            # 두 번째 OH 흡착 에너지
            OH_Ocovered.append(energies_o_ohcus[i] + dgoh - energies_v_o[i] - oh_gas_energy)

        return ueff_values, O_Ocovered, OH_Ocovered
    else:
        return None, None, None


def plot_o_minus_oh_vs_ueff(ueff_values, O_Ocovered, OH_Ocovered):
    # 차이값 계산: Oads - OHads
    D_Ocovered = [o - oh for o, oh in zip(O_Ocovered, OH_Ocovered)]

    plt.figure(figsize=(4, 4))

    plt.plot(ueff_values, D_Ocovered, color='black', marker='o', label=r'$\Delta G_{O} - \Delta G_{OH}$')

    plt.xlabel(r'$U_{\mathrm{eff}}$ (eV)', fontsize=12)
    plt.ylabel(r'$\Delta G_{O+O} - \Delta G_{O+OH}$ (eV)', fontsize=12)
    plt.xlim(-0.5, 4.5)
    
    # 데이터 포인트에 값 표시
    for i, (x, y) in enumerate(zip(ueff_values, D_Ocovered)):
        if i == 0 or i == 2 or i == 1:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", xytext=(0,-14), ha='center', va='center')
        else:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", xytext=(0,10), ha='center', va='center')

    plt.ylim(0.8, 1.8)
    plt.axhline(1.327, color='blue', linewidth=1.0, linestyle='--')
    plt.axhline(1.573, color='blue', linewidth=1.0, linestyle='--')
    plt.scatter(2.5, 1.327, color='blue', marker='s')
    plt.scatter(3.9, 1.573, color='blue', marker='s')
    # plt.axhspan(1.23, 1.60, color='lightsteelblue', zorder=0, alpha=0.5)    # mistyrose, gainsboro

    plt.tight_layout()
    plt.savefig('OminusOH_adsorption_energies_vs_Ueff.png', dpi=300, bbox_inches='tight')
    plt.show()

    # 결과 저장 (차이값)
    with open('OminusOH_adsorption_energies_results.txt', 'w') as f:
        f.write("Ueff (eV)\t$\Delta G_{O} - \Delta G_{OH}$ (top, O-covered) (eV)\n")
        f.write("-" * 60 + "\n")
        for i, u in enumerate(ueff_values):
            f.write(f"{u}\t\t{D_Ocovered[i]:.6f}\n")


if __name__ == "__main__":
    result = calculate_o_and_oh_adsorption_energies()
    if result[0] is None:
        print("필요한 에너지를 모두 불러오지 못했습니다. 경로와 파일을 확인하세요.")
    else:
        ueff_values, O_Ocovered, OH_Ocovered = result
        plot_o_minus_oh_vs_ueff(ueff_values, O_Ocovered, OH_Ocovered)


