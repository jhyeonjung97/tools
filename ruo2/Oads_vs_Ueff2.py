from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

root = "~/Desktop/3_RuO2/1_RuO2_Ueff"

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
    v_o_path = root_path / "3_O_V"
    o_o_path = root_path / "5_O_O"
    o_ohcus_path = root_path / "4_O_OH"
    
    ueff_values = [0, 1, 2, 3, 4]
    energies_v_o = []
    energies_o_o = []
    energies_o_ohcus = []
    
    # V_O 에너지 추출
    for ueff in ueff_values:
        json_path = v_o_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_v_o.append(energy)
    
    # O_O 에너지 추출
    for ueff in ueff_values:
        json_path = o_o_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_o_o.append(energy)
    
    # O_OHcus 에너지 추출
    for ueff in ueff_values:
        json_path = o_ohcus_path / f"{ueff}_" / "final_with_calculator.traj"
        energy = get_energy_from_json(json_path)
        if energy is not None:
            energies_o_ohcus.append(energy)
    
    if len(energies_v_o) == len(energies_o_o) == len(energies_o_ohcus) == len(ueff_values):
        o_adsorption_energies_2 = []  # O_O - V_O (두 번째 O 흡착)
        oh_adsorption_energies_2 = []  # O_OHcus - V_O (두 번째 OH 흡착)
        
        for i, ueff in enumerate(ueff_values):
            # 두 번째 O 흡착 에너지: O가 이미 흡착된 표면에서 추가 O 흡착
            o_gas_energy = gh2o - gh2
            o_ads_energy_2 = energies_o_o[i] + dgo - energies_v_o[i] - o_gas_energy
            o_adsorption_energies_2.append(o_ads_energy_2)
            
            # 두 번째 OH 흡착 에너지: O가 이미 흡착된 표면에서 추가 OH 흡착
            oh_gas_energy = gh2o - gh2/2
            oh_ads_energy_2 = energies_o_ohcus[i] + dgoh - energies_v_o[i] - oh_gas_energy
            oh_adsorption_energies_2.append(oh_ads_energy_2)
        
        return ueff_values, o_adsorption_energies_2, oh_adsorption_energies_2
    else:
        return None, None, None

def plot_adsorption_energies_vs_ueff(ueff_values, o_energies_2, oh_energies_2):    
    plt.figure(figsize=(4, 4))
    
    # 두 번째 O 흡착 에너지: O-covered surface에서 O 흡착
    plt.plot(ueff_values, o_energies_2, color='black', marker='o', markeredgecolor='black', markerfacecolor='black', label=r'$\Delta G_{O}$')
    
    # 두 번째 OH 흡착 에너지: O-covered surface에서 OH 흡착
    plt.plot(ueff_values, oh_energies_2, color='black', marker='o', markeredgecolor='black', markerfacecolor='white', label=r'$\Delta G_{OH}$')
    
    plt.xlabel(r'$U_{\mathrm{eff}}$ (eV)')
    plt.ylabel(r'$\Delta G_{O}$ or $\Delta G_{OH}$ (eV)')
    plt.xlim(-0.5, 4.5)
    plt.legend()
    
    # 데이터 포인트에 값 표시
    for i, (x, y) in enumerate(zip(ueff_values, o_energies_2)):
        plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", xytext=(0,-13), ha='center', va='center')
    
    for i, (x, y) in enumerate(zip(ueff_values, oh_energies_2)):
        plt.annotate(f'{y:.2f}', (x, y), textcoords="offset points", xytext=(0,10), ha='center', va='center')
    
    # y축 범위를 두 데이터 모두 포함하도록 설정
    all_values = o_energies_2 + oh_energies_2
    plt.ylim(0.5, 2.5)
    plt.yticks(np.arange(0.5, 2.6, 0.4))

    plt.tight_layout()
    plt.savefig('O_OH_adsorption_energies_vs_Ueff.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 결과를 텍스트 파일로 저장
    with open('O_OH_adsorption_energies_results.txt', 'w') as f:
        f.write("Ueff (eV)\t$\Delta G_{O}$ (top, O-covered) (eV)\t$\Delta G_{OH}$ (top, O-covered) (eV)\n")
        f.write("-" * 80 + "\n")
        for ueff, o_energy_2, oh_energy_2 in zip(ueff_values, o_energies_2, oh_energies_2):
            f.write(f"{ueff}\t\t{o_energy_2:.6f}\t\t\t{oh_energy_2:.6f}\n")

if __name__ == "__main__":
    ueff_values, o_energies_2, oh_energies_2 = calculate_adsorption_energies()
    
    if ueff_values is not None and o_energies_2 is not None and oh_energies_2 is not None:
        plot_adsorption_energies_vs_ueff(ueff_values, o_energies_2, oh_energies_2)
