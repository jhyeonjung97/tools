from ase.io import read
from pathlib import Path

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

calmol = 23.061
ru = -9.243204045
ru2 = ru + 0.8*2
ru3 = ru2 + 0.249
# ruo2h2 = ru + 0.68*4 + 2*gh2o - gh2
# # ruo2h2 = ru + 0.8*2 + 0.249 + 0.86 + 2*gh2o - gh2\

# ruo2 = ru + 0.8*2 + 0.249 + 0.86 + 2*gh2o - 2*gh2 + 1.5*2
# ruo4 = ru + 4*gh2o - 4*h2 + 1.04*8
# ruo4_ = ru + 0.8*2 + 0.249 + 0.86 + 2*gh2o - 2*gh2 + 1.5*2 + 1.6 + 0.99

# ru2 = ru + 21.000/calmol
# ruo4 = ru + 0.8*2 + 0.249 + 0.86 + 2*gh2o - 2*gh2 + 1.5*2 + 1.6
# print(ru + 0.8*2 + 0.249 + 0.86 + 2*gh2o - 2*gh2 + 1.5*2 + 1.6 + 0.99)
# print(ru + 4*gh2o - 4*h2 + 48.000/calmol)
ruo4_ = ru + 4*gh2o - 4*gh2+ 48.000/calmol
ruo4__ = ru + 4*gh2o - 4*gh2 + 61.600/calmol
h2ruo5 = ru + 5*gh2o - 4*gh2 + 81.600/calmol
# hruo5 = ru + 5*go + gh2/2 - 66.300/calmol

def get_energy_from_json(json_path):
    """ASE를 사용하여 JSON 파일에서 에너지를 추출합니다."""
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def read_dissolution_energies():
    """Dissolution 폴더에서 에너지 데이터를 읽습니다."""
    dissolution_root = "/Users/jiuy97/Desktop/3_RuO2/7_Dissolution_110"
    dissolution_path = Path(dissolution_root)
    
    # 메인 폴더들
    main_folders = ["0_RuO2", "1_ReRuO2_Top1", "2_ReRuO2_Top2", "3_ReRuO2_Brg1", "4_ReRuO2_Brg2"]
    # 하위 폴더들
    sub_folders = ["1_Ru_Top", "2_RuO2_Top", "3_Ru_Brg", "4_RuO2_Brg", "5_O_Brg", "6_O2_Brg"]
    
    dissolution_data = {}
    
    for main_folder in main_folders:
        dissolution_data[main_folder] = {}
        for sub_folder in sub_folders:
            json_path = dissolution_path / main_folder / sub_folder / "final_with_calculator.json"
            energy = get_energy_from_json(json_path)
            dissolution_data[main_folder][sub_folder] = energy
            print(f"Dissolution - {main_folder}/{sub_folder}: {energy}")
    
    return dissolution_data

def read_intermediates_energies():
    """Intermediates 폴더에서 에너지 데이터를 읽습니다."""
    intermediates_root = "/Users/jiuy97/Desktop/3_RuO2/6_Intermediates_110"
    intermediates_path = Path(intermediates_root)
    
    # 메인 폴더들
    main_folders = ["0_RuO2", "1_ReRuO2_Top1", "2_ReRuO2_Top2", "3_ReRuO2_Brg1", "4_ReRuO2_Brg2"]
    # 하위 폴더 (Intermediates에서는 6_O_O만 있음)
    sub_folder = "6_O_O"
    
    intermediates_data = {}
    
    for main_folder in main_folders:
        json_path = intermediates_path / main_folder / sub_folder / "final_with_calculator.json"
        energy = get_energy_from_json(json_path)
        intermediates_data[main_folder] = energy
        print(f"Intermediates - {main_folder}/{sub_folder}: {energy}")
    
    return intermediates_data

def print_energies_summary(dissolution_data, intermediates_data):
    """에너지 데이터를 화면에 출력합니다."""
    print("=== Energy Summary ===")
    print()
    
    # Dissolution 데이터 출력
    print("Dissolution Energies:")
    print("-" * 50)
    for main_folder, energies in dissolution_data.items():
        print(f"\n{main_folder}:")
        for sub_folder, energy in energies.items():
            if energy is not None:
                print(f"  {sub_folder}: {energy:.6f} eV")
            else:
                print(f"  {sub_folder}: Not found")
    
    print("\n" + "="*50)
    
    # Intermediates 데이터 출력
    print("\nIntermediates Energies:")
    print("-" * 50)
    for main_folder, energy in intermediates_data.items():
        if energy is not None:
            print(f"{main_folder}/6_O_O: {energy:.6f} eV")
        else:
            print(f"{main_folder}/6_O_O: Not found")

def calculate_energy_differences(dissolution_data, intermediates_data):
    """각 표면에 대해 에너지 차이를 계산합니다."""
    print("=== Energy Differences Analysis ===")
    print()
    
    main_folders = ["0_RuO2", "1_ReRuO2_Top1", "2_ReRuO2_Top2", "3_ReRuO2_Brg1", "4_ReRuO2_Brg2"]
    
    results = {}
    
    for main_folder in main_folders:
        print(f"--- {main_folder} ---")
        
        # 4_RuO2_Brg 에너지 (기준)
        ruo2_brg_energy = dissolution_data[main_folder]["4_RuO2_Brg"]
        
        if ruo2_brg_energy is None:
            print(f"Warning: 4_RuO2_Brg energy not found for {main_folder}")
            continue
            
        results[main_folder] = {}
        
        # (1) 6_O_O (Intermediates) → 4_RuO2_Brg
        o_o_energy = intermediates_data[main_folder]
        if o_o_energy is not None:
            # diff_1 = -(o_o_energy + 3*gh2o - ruo2_brg_energy - h2ruo5 - 2*gh2)/4
            # diff_1 = -(o_o_energy + 2*gh2o- ruo2_brg_energy - 2*gh2 - ruo4_)/3
            diff_1 = -(o_o_energy + 2*gh2o- ruo2_brg_energy - 2*gh2 - ruo4__)/3
            # diff_1 = -(o_o_energy - ruo2_brg_energy - ruo2)/2
            # diff_1 = -(o_o_energy + gh2 - ruo2_brg_energy - ruo2h2)
            results[main_folder]["6_O_O_to_4_RuO2_Brg"] = diff_1
            print(f"6_O_O → 4_RuO2_Brg: {diff_1:.6f} eV")
        else:
            print(f"Warning: 6_O_O energy not found for {main_folder}")
            results[main_folder]["6_O_O_to_4_RuO2_Brg"] = None
        
        # (2) 5_O_Brg → 4_RuO2_Brg
        o_brg_energy = dissolution_data[main_folder]["5_O_Brg"]
        if o_brg_energy is not None:
            diff_2 = -(o_brg_energy - ruo2_brg_energy - ru3)
            # diff_2 = -(o_brg_energy - ruo2_brg_energy - ru3)
            # diff_2 = -(o_brg_energy + gh2 - ruo2_brg_energy - ru2 -gh2o)
            results[main_folder]["5_O_Brg_to_4_RuO2_Brg"] = diff_2
            print(f"5_O_Brg → 4_RuO2_Brg: {diff_2:.6f} eV")
        else:
            print(f"Warning: 5_O_Brg energy not found for {main_folder}")
            results[main_folder]["5_O_Brg_to_4_RuO2_Brg"] = None
        
        # (3) 6_O2_Brg → 4_RuO2_Brg
        o2_brg_energy = dissolution_data[main_folder]["6_O2_Brg"]
        if o2_brg_energy is not None:
            # diff_3 = -(o2_brg_energy - ruo2_brg_energy - ru3)/3
            diff_3 = -(o2_brg_energy - ruo2_brg_energy - ru)
            results[main_folder]["6_O2_Brg_to_4_RuO2_Brg"] = diff_3
            print(f"6_O2_Brg → 4_RuO2_Brg: {diff_3:.6f} eV")
        else:
            print(f"Warning: 6_O2_Brg energy not found for {main_folder}")
            results[main_folder]["6_O2_Brg_to_4_RuO2_Brg"] = None
        
        print()
    
    return results

def print_energy_differences_summary(energy_differences):
    """에너지 차이 결과를 화면에 출력합니다."""
    print("\n=== Energy Differences Summary ===")
    print()
    print(f"{'Surface':<20} {'6_O_O→4_RuO2_Brg':<20} {'5_O_Brg→4_RuO2_Brg':<20} {'6_O2_Brg→4_RuO2_Brg':<20}")
    print("-" * 80)
    
    for surface, differences in energy_differences.items():
        diff1 = differences.get("6_O_O_to_4_RuO2_Brg", "N/A")
        diff2 = differences.get("5_O_Brg_to_4_RuO2_Brg", "N/A")
        diff3 = differences.get("6_O2_Brg_to_4_RuO2_Brg", "N/A")
        
        if diff1 != "N/A":
            diff1 = f"{diff1:.6f}"
        if diff2 != "N/A":
            diff2 = f"{diff2:.6f}"
        if diff3 != "N/A":
            diff3 = f"{diff3:.6f}"
            
        print(f"{surface:<20} {diff1:<20} {diff2:<20} {diff3:<20}")

def main():
    """메인 함수"""
    print("=== RuO2 Dissolution Energy Analysis ===")
    print()
    
    print("Reading dissolution energies...")
    dissolution_data = read_dissolution_energies()
    print()
    
    print("Reading intermediates energies...")
    intermediates_data = read_intermediates_energies()
    print()
    
    print_energies_summary(dissolution_data, intermediates_data)
    print()
    
    print("Calculating energy differences...")
    energy_differences = calculate_energy_differences(dissolution_data, intermediates_data)
    
    print_energy_differences_summary(energy_differences)
    print()
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()
