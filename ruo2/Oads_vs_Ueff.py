from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path

root = "~/Desktop/3_RuO2/1_RuO2_110"

# gas
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

def calculate_adsorption_energies():    
    root_path = Path(root).expanduser()
    v_v_path = root_path / "1_V_V"
    v_o_path = root_path / "3_V_O"
    o_o_path = root_path / "6_O_O"
    vo_o_path = root_path / "7_VO_O"
    
    ueff_values = [0, 1, 2, 3, 4]
    energies_v_v = []
    energies_v_o = []
    energies_o_o = []
    energies_vo_o = []
        
    # V_V 에너지 추출
    for ueff in ueff_values:
        json_path = v_v_path / f"{ueff}_" / "final_with_calculator.json"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_v_v.append(energy)
    
    # V_O 에너지 추출
    for ueff in ueff_values:
        json_path = v_o_path / f"{ueff}_" / "final_with_calculator.json"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_v_o.append(energy)
    
    # O_O 에너지 추출
    for ueff in ueff_values:
        json_path = o_o_path / f"{ueff}_" / "final_with_calculator.json"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_o_o.append(energy)
    
    # VO_O 에너지 추출
    for ueff in ueff_values:
        json_path = vo_o_path / f"{ueff}_" / "final_with_calculator.json"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_vo_o.append(energy)
    
    if len(energies_v_v) == len(energies_v_o) == len(energies_o_o) == len(energies_vo_o) == len(ueff_values):
        o_adsorption_energies_1 = []  # V_O - V_V
        o_adsorption_energies_2 = []  # O_O - V_O
        o_adsorption_energies_3 = []  # VO_O - O_O
        
        for i, ueff in enumerate(ueff_values):
            # 첫 번째 O 흡착 에너지: clean surface에서 O 흡착
            o_gas_energy = gh2o - gh2
            o_ads_energy_1 = energies_v_o[i] + dgo - energies_v_v[i] - o_gas_energy
            o_adsorption_energies_1.append(o_ads_energy_1)
            print(f"Ueff={ueff}: dG_O (clean) = {o_ads_energy_1:.6f} eV")
            
            # 두 번째 O 흡착 에너지: O가 이미 흡착된 표면에서 추가 O 흡착
            o_ads_energy_2 = energies_o_o[i] + dgo - energies_v_o[i] - o_gas_energy
            o_adsorption_energies_2.append(o_ads_energy_2)
            print(f"Ueff={ueff}: dG_O (O-covered) = {o_ads_energy_2:.6f} eV")
            
            # 세 번째 O 흡착 에너지: VO_O에서 O_O로의 O 흡착
            o_ads_energy_3 = energies_o_o[i] + dgo - energies_vo_o[i] - o_gas_energy
            o_adsorption_energies_3.append(o_ads_energy_3)
            print(f"Ueff={ueff}: dG_O (VO_O) = {o_ads_energy_3:.6f} eV")
        
        return ueff_values, o_adsorption_energies_1, o_adsorption_energies_2, o_adsorption_energies_3
    else:
        return None, None, None, None

def plot_o_adsorption_energies_vs_ueff(ueff_values, o_energies_1, o_energies_2, o_energies_3):    
    plt.figure(figsize=(6, 5))
    
    # 첫 번째 O 흡착 에너지: clean surface에서 O 흡착
    plt.plot(ueff_values, o_energies_1, 'bo-', linewidth=2, markersize=8, label=r'$\Delta G_{O}$ (top)')
    
    # 두 번째 O 흡착 에너지: O-covered surface에서 O 흡착
    plt.plot(ueff_values, o_energies_2, 'rs-', linewidth=2, markersize=8, label=r'$\Delta G_{O}$ (top, O-covered)')
    
    # 세 번째 O 흡착 에너지: VO_O에서 O_O로의 O 흡착
    plt.plot(ueff_values, o_energies_3, 'g^-', linewidth=2, markersize=8, label=r'$\Delta G_{O}$ (brg, O-covered)')
    
    plt.xlabel(r'$U_{\mathrm{eff}}$ (eV)', fontsize=14)
    plt.ylabel(r'$\Delta G_{O}$ (eV)', fontsize=14)
    plt.xlim(-0.5, 4.5)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    
    # 데이터 포인트에 값 표시
    for i, (x, y) in enumerate(zip(ueff_values, o_energies_1)):
        plt.annotate(f'{y:.3f}', (x, y), textcoords="offset points", 
                    xytext=(0,-15), ha='center', fontsize=9)
    
    for i, (x, y) in enumerate(zip(ueff_values, o_energies_2)):
        plt.annotate(f'{y:.3f}', (x, y), textcoords="offset points", 
                    xytext=(0,10), ha='center', fontsize=9)

    for i, (x, y) in enumerate(zip(ueff_values, o_energies_3)):
        plt.annotate(f'{y:.3f}', (x, y), textcoords="offset points", 
                    xytext=(0,10), ha='center', fontsize=9)
    
    # plt.ylim(-1.5, 1.0)
    
    plt.tight_layout()
    plt.savefig('O_adsorption_energies_vs_Ueff.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 결과를 텍스트 파일로 저장
    with open('O_adsorption_energies_results.txt', 'w') as f:
        f.write("Ueff (eV)\t$\Delta G_{O}$ (top) (eV)\t$\Delta G_{O}$ (top, O-covered) (eV)\t$\Delta G_{O}$ (brg, O-covered) (eV)\n")
        f.write("-" * 100 + "\n")
        for ueff, o_energy_1, o_energy_2, o_energy_3 in zip(ueff_values, o_energies_1, o_energies_2, o_energies_3):
            f.write(f"{ueff}\t\t{o_energy_1:.6f}\t\t\t{o_energy_2:.6f}\t\t\t{o_energy_3:.6f}\n")

if __name__ == "__main__":
    ueff_values, o_energies_1, o_energies_2, o_energies_3 = calculate_adsorption_energies()
    
    if ueff_values is not None and o_energies_1 is not None and o_energies_2 is not None and o_energies_3 is not None:
        plot_o_adsorption_energies_vs_ueff(ueff_values, o_energies_1, o_energies_2, o_energies_3)
