from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path

# 사용자 환경에 맞춰 경로를 변경하세요
root = "~/Desktop/3_RuO2/1_RuO2_110"

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
    """O와 OH 흡착 자유에너지를 Ueff에 대해 계산.

    반환:
        ueff_values, (O_clean, O_Ocovered, O_from_VO), (OH_clean, OH_Ocovered, OH_from_VO)
    """
    root_path = Path(root).expanduser()

    # 경로 매핑 (기존 두 스크립트와 동일한 폴더 구조 사용)
    v_v_path = root_path / "1_V_V"
    v_oh_path = root_path / "2_V_OH"
    v_o_path = root_path / "3_V_O"
    o_ohcus_path = root_path / "4_O_OH"
    o_o_path = root_path / "5_O_O"
    vo_o_path = root_path / "6_VO_O"
    ho_o_path = root_path / "7_HO_O"

    ueff_values = [0, 1, 2, 3, 4]

    energies_v_v = []
    energies_v_oh = []
    energies_v_o = []
    energies_o_ohcus = []
    energies_o_o = []
    energies_vo_o = []
    energies_ho_o = []

    # 에너지 로드 유틸
    def load_series(base_path, label):
        series = []
        for u in ueff_values:
            json_path = base_path / f"{u}_" / "final_with_calculator.json"
            e = get_energy_from_json(json_path)
            if e is None:
                print(f"Skip {label} Ueff={u} (missing/failed)")
                return None
            series.append(e)
        return series

    energies_v_v = load_series(v_v_path, "V_V")
    energies_v_oh = load_series(v_oh_path, "V_OH")
    energies_v_o = load_series(v_o_path, "V_O")
    energies_o_ohcus = load_series(o_ohcus_path, "O_OHcus")
    energies_o_o = load_series(o_o_path, "O_O")
    energies_vo_o = load_series(vo_o_path, "VO_O")
    energies_ho_o = load_series(ho_o_path, "HO_O")

    if any(x is None for x in [energies_v_v, energies_v_oh, energies_v_o, energies_o_ohcus, energies_o_o, energies_vo_o, energies_ho_o]):
        return None

    # O 흡착 에너지
    O_clean = []       # V_O - V_V
    O_Ocovered = []    # O_O - V_O
    O_from_VO = []     # O_O - VO_O

    # OH 흡착 에너지
    OH_clean = []      # V_OH - V_V
    OH_Ocovered = []   # O_OHcus - V_O
    OH_from_VO = []    # HO_O - VO_O

    for i, ueff in enumerate(ueff_values):
        # 기준 가스 상
        o_gas_energy = gh2o - gh2
        oh_gas_energy = gh2o - gh2/2

        # O
        O_clean.append(energies_v_o[i] + dgo - energies_v_v[i] - o_gas_energy)
        O_Ocovered.append(energies_o_o[i] + dgo - energies_v_o[i] - o_gas_energy)
        O_from_VO.append(energies_o_o[i] + dgo - energies_vo_o[i] - o_gas_energy)

        # OH
        OH_clean.append(energies_v_oh[i] + dgoh - energies_v_v[i] - oh_gas_energy)
        OH_Ocovered.append(energies_o_ohcus[i] + dgoh - energies_v_o[i] - oh_gas_energy)
        OH_from_VO.append(energies_ho_o[i] + dgoh - energies_vo_o[i] - oh_gas_energy)

    return (
        ueff_values,
        (O_clean, O_Ocovered, O_from_VO),
        (OH_clean, OH_Ocovered, OH_from_VO),
    )


def plot_o_minus_oh_vs_ueff(ueff_values, O_series, OH_series):
    O_clean, O_Ocovered, O_from_VO = O_series
    OH_clean, OH_Ocovered, OH_from_VO = OH_series

    # 차이값 계산: Oads - OHads
    D_clean = [o - oh for o, oh in zip(O_clean, OH_clean)]
    D_Ocovered = [o - oh for o, oh in zip(O_Ocovered, OH_Ocovered)]
    D_from_VO = [o - oh for o, oh in zip(O_from_VO, OH_from_VO)]

    plt.figure(figsize=(5.5, 5))

    # # 수평선 추가
    # plt.axhline(y=1.22, color='black', linestyle='-', alpha=0.7, linewidth=1)
    # plt.axhline(y=1.31, color='black', linestyle='--', alpha=0.7, linewidth=1)

    # plt.text(0.3, 1.22+0.01, r'sRuO$_2$', fontsize=10, color='black')
    # plt.text(0.3, 1.31+0.01, r'cRuO$_2$', fontsize=10, color='black')

    plt.plot(ueff_values, D_clean, 'bo-', linewidth=2, markersize=8, label=r'$\Delta G_{O} - \Delta G_{OH}$ (top, clean)')
    plt.plot(ueff_values, D_Ocovered, 'rs-', linewidth=2, markersize=8, label=r'$\Delta G_{O} - \Delta G_{OH}$ (top, O-covered)')
    plt.plot(ueff_values, D_from_VO, 'g^-', linewidth=2, markersize=8, label=r'$\Delta G_{O} - \Delta G_{OH}$ (brg, O-covered)')

    plt.xlabel(r'$U_{\mathrm{eff}}$ (eV)', fontsize=14)
    plt.ylabel(r'$\Delta G_{O} - \Delta G_{OH}$ (eV)', fontsize=14)
    plt.xlim(-0.5, 4.5)
    plt.ylim(0.8, 1.7)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # 데이터 포인트에 값 표시
    for i, (x, y) in enumerate(zip(ueff_values, D_clean)):
        if i != 2 and i != 3:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,-15), ha='center', fontsize=9)
        else:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,10), ha='center', fontsize=9)
    
    for i, (x, y) in enumerate(zip(ueff_values, D_Ocovered)):
        if i != 1:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,10), ha='center', fontsize=9)
        else:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,-15), ha='center', fontsize=9)

    for i, (x, y) in enumerate(zip(ueff_values, D_from_VO)):
        if i != 2 and i != 3:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,10), ha='center', fontsize=9)
        else:
            plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", 
                        xytext=(0,-15), ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig('OminusOH_adsorption_energies_vs_Ueff.png', dpi=300, bbox_inches='tight')
    plt.show()

    # 결과 저장 (차이값)
    with open('OminusOH_adsorption_energies_results.txt', 'w') as f:
        f.write("Ueff (eV)\t(D_O-D_OH)_clean\t(D_O-D_OH)_Ocov\t(D_O-D_OH)_fromVO\n")
        f.write("-" * 80 + "\n")
        for i, u in enumerate(ueff_values):
            f.write(f"{u}\t{D_clean[i]:.6f}\t{D_Ocovered[i]:.6f}\t{D_from_VO[i]:.6f}\n")


if __name__ == "__main__":
    result = calculate_o_and_oh_adsorption_energies()
    if result is None:
        print("필요한 에너지를 모두 불러오지 못했습니다. 경로와 파일을 확인하세요.")
    else:
        ueff_values, O_series, OH_series = result
        plot_o_minus_oh_vs_ueff(ueff_values, O_series, OH_series)


