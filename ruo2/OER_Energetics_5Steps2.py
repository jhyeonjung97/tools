from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# 폰트 설정: Arial
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# 사용자 환경에 맞춰 경로를 변경하세요
root = "~/Desktop/3_RuO2"

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
gh = gh2/2
go = gh2o - gh2
goh = gh2o - gh
goo = 2*gh2o - 2*gh2
gooh = 2*gh2o - 3*gh

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
dgoo = dgooh - dgh

def get_energy_from_json(json_path):
    """JSON 파일에서 에너지를 읽어오는 함수"""
    try:
        atoms = read(json_path)
        energy = atoms.get_total_energy()
        return energy
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return None

def calculate_oer_energies():
    """14가지 OER 경로의 반응 에너지를 계산하고 표 형태로 정리하는 함수"""
    root_path = Path(root).expanduser()

    gibbs_correction_o_v = 0.058092 - 0.033706
    gibbs_correction_o_oh = 0.439910 - 0.128005
    gibbs_correction_o_o = 0.156862 - 0.100681
    gibbs_correction_oh_oo = 0.502349 - 0.144209
    gibbs_correction_v_oh = 0.364538 - 0.074708

    # 표 형태로 데이터를 저장할 리스트
    oer_data = []
    all_paths = {}

    try:
        oer_path = root_path / "2_RuO2_OER" / "7_RuO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_o_ooh = get_energy_from_json(oer_path / "4_O_OOH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_o_ooh, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_o_ooh - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2(a)',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'O_OOH',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path1'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 1: Failed to load some energies")
    except Exception as e:
        print(f"Path 1: Error - {e}")

    try:
        oer_path = root_path / "2_RuO2_OER" / "7_RuO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh
        
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5

            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2(d)',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path2'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 2: Failed to load some energies")
    except Exception as e:
        print(f"Path 2: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "3_ReRuO2"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path3'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 3: Failed to load some energies")
    except Exception as e:
        print(f"Path 3: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "4_RuO2_1%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2_1%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path4'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 4: Failed to load some energies")
    except Exception as e:
        print(f"Path 4: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "5_RuO2_2%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh
        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'RuO2_2%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path5'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 5: Failed to load some energies")
    except Exception as e:
        print(f"Path 5: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "6_ReRuO2_1%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2_1%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path6'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 6: Failed to load some energies")
    except Exception as e:
        print(f"Path 6: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "7_ReRuO2_2%"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2_2%',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'OH_OO',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path7'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 7: Failed to load some energies")
    except Exception as e:
        print(f"Path 7: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "9_double_cell"
        energy_o_v = get_energy_from_json(oer_path / "0_" / "0_" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "1_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "2_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "3_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "4_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            # step4 = (energy_v_oh) - (energy_oh_oo - goo) + 4.92
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2-V1',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'O_OOH',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path8'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 8: Failed to load some energies")
    except Exception as e:
        print(f"Path 8: Error - {e}")

    try:
        oer_path = root_path / "3_ReRuO2_OER" / "9_double_cell"
        energy_o_v = get_energy_from_json(oer_path / "0_" / "1_" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "5_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "6_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "7_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "8_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            print('d2')
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            # step4 = (energy_v_oh) - (energy_oh_oo - goo) + 4.92
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2-V2',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'O_OOH',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path9'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 9: Failed to load some energies")
    except Exception as e:
        print(f"Path 9: Error - {e}")

    try:
        oer_path = root_path / "5_ReRuO2_alloy" / "3_OER"
        energy_o_v = get_energy_from_json(oer_path / "1_O_V" / "final_with_calculator.traj") + gibbs_correction_o_v
        energy_o_oh = get_energy_from_json(oer_path / "2_O_OH" / "final_with_calculator.traj") + gibbs_correction_o_oh
        energy_o_o = get_energy_from_json(oer_path / "3_O_O" / "final_with_calculator.traj") + gibbs_correction_o_o
        energy_oh_oo = get_energy_from_json(oer_path / "5_OO_OH" / "final_with_calculator.traj") + gibbs_correction_oh_oo
        energy_v_oh = get_energy_from_json(oer_path / "6_V_OH" / "final_with_calculator.traj") + gibbs_correction_v_oh

        if all(e is not None for e in [energy_o_v, energy_o_oh, energy_o_o, energy_oh_oo, energy_v_oh]):
            step1 = (energy_o_oh - goh) - (energy_o_v)
            step2 = (energy_o_o - go) - (energy_o_oh - goh)
            step3 = (energy_oh_oo - goh) - (energy_o_o)
            # step4 = (energy_v_oh) - (energy_oh_oo - goo) + 4.92
            step5 = (energy_o_v - go) - (energy_v_oh - goh)
            step4 = 4.92 - step1 - step2 - step3 - step5
            
            # 표 형태로 데이터 추가
            oer_data.append({
                'surface': 'ReRuO2-alloy',
                'int1': 'O_V',
                'int2': 'O_OH', 
                'int3': 'O_O',
                'int4': 'O_OOH',
                'int5': 'V_OH',
                'step1': step1,
                'step2': step2,
                'step3': step3,
                'step4': step4,
                'step5': step5
            })
            all_paths['Path10'] = [step1, step2, step3, step4, step5]
            print(f"Step1: {step1:.3f} eV, Step2: {step2:.3f} eV, Step3: {step3:.3f} eV, Step4: {step4:.3f} eV, Step5: {step5:.3f} eV")
        else:
            print("Path 10: Failed to load some energies")
    except Exception as e:
        print(f"Path 10: Error - {e}")

    return all_paths, oer_data

def print_oer_table(oer_data):
    """OER 데이터를 표 형태로 출력하는 함수"""
    if not oer_data:
        print("No OER data to display")
        return
    
    print("\n" + "="*120)
    print("OER Energetics Summary Table")
    print("="*120)
    
    # 헤더 출력
    header = f"{'Surface':<12} {'Int1':<8} {'Int2':<8} {'Int3':<8} {'Int4':<8} {'Int5':<8} {'Step1':<10} {'Step2':<10} {'Step3':<10} {'Step4':<10} {'Step5':<10} {'Max':<10} {'Overpotential':<12}"
    print(header)
    print("-"*120)
    
    # 데이터 출력
    for i, data in enumerate(oer_data, 1):
        max_energy = max(data['step1'], data['step2'], data['step3'], data['step4'], data['step5'])
        overpotential = max_energy - 1.23
        
        row = f"{data['surface']:<12} {data['int1']:<8} {data['int2']:<8} {data['int3']:<8} {data['int4']:<8} {data['int5']:<8} " \
              f"{data['step1']:<10.3f} {data['step2']:<10.3f} {data['step3']:<10.3f} {data['step4']:<10.3f} {data['step5']:<10.3f} " \
              f"{max_energy:<10.3f} {overpotential:<12.3f}"
        print(row)
    
    print("-"*120)
    
    # 최적 경로 찾기
    best_path = min(oer_data, key=lambda x: max(x['step1'], x['step2'], x['step3'], x['step4'], x['step5']))
    best_max = max(best_path['step1'], best_path['step2'], best_path['step3'], best_path['step4'], best_path['step5'])
    best_overpotential = best_max - 1.23
    
    print(f"\nBest Path: {best_path['surface']} - {best_path['int1']} → {best_path['int2']} → {best_path['int3']} → {best_path['int4']} → {best_path['int5']}")
    print(f"Maximum Step Energy: {best_max:.3f} eV")
    print(f"Overpotential: {best_overpotential:.3f} eV")
    print("="*120)

def plot_individual_energetics_diagrams(all_paths, oer_data):
    """각 path에 대해 개별 energetics diagram을 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting individual diagrams")
        return
    
    # 각 path별로 개별 그래프 생성
    for i, (path_name, energies) in enumerate(all_paths.items()):
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))

        # 데이터에서 해당 path 정보 찾기 (인덱스로 찾기)
        path_data = oer_data[i] if i < len(oer_data) else None
        if not path_data:
            continue
            
        # x축: 반응 단계 (1, 2, 3, 4, 5)
        steps = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]

        energy0 = 0.0
        energy1 = path_data['step1']
        energy2 = path_data['step2'] + energy1
        energy3 = path_data['step3'] + energy2
        energy4 = path_data['step4'] + energy3
        energy5 = path_data['step5'] + energy4
        surface = path_data['surface']
        
        # y축: step energies
        step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4, energy5, energy5]
        
        # 선 그래프로 그리기
        ax.plot(steps, step_energies, '-', color='black', linewidth=2)

        # 그래프 설정
        ax.set_xticks([])
        ax.set_xlabel('Reaction Coordinate', fontsize=12)
        ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
        
        # y축 범위 설정 (최소값과 최대값에 여유 공간 추가)
        y_min = min(step_energies) - 0.2
        y_max = max(step_energies) + 0.3
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(0, 6)
                
        # 최대 에너지와 과전위 정보 추가
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'], path_data['step5'])
        overpotential = max_energy - 1.23
        ax.text(0.02, 0.97, f'Surface: {surface}\nΔG1: {path_data["step1"]:.2f} eV\nΔG2: {path_data["step2"]:.2f} eV\nΔG3: {path_data["step3"]:.2f} eV\nΔG4: {path_data["step4"]:.2f} eV\nΔG5: {path_data["step5"]:.2f} eV\nOverpotential: {overpotential:.2f} V', 
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        # 파일명에 path 정보 포함하여 저장
        filename = f'energetics_diagram_5steps_{path_name.lower()}_{path_data["surface"].lower()}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        # plt.show()
        plt.close()
        
        print(f"Saved: {filename}")

def plot_all_pathways_combined_RuO2(all_paths, oer_data):
    """RuO2, ReRuO2, RuO2_1%, RuO2_2% 케이스를 하나의 그래프에 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting combined pathways")
        return
    
    # 필터링할 surface 목록
    target_surfaces = ['RuO2(d)', 'ReRuO2', 'RuO2_1%', 'RuO2_2%']
    names = ['RuO$_2$', 'Re(brg)-RuO$_2$', 'RuO$_2$(+1%)', 'RuO$_2$(+2%)'] #'Re(brg)-RuO$_2$(+1%)', 'Re(brg)-RuO$_2$(+2%)']

    # 필터링된 데이터만 추출
    filtered_data = []
    
    for path_data in oer_data:
        if path_data['surface'] in target_surfaces:
            filtered_data.append(path_data)
    
    if not filtered_data:
        print("No matching data found for RuO2, ReRuO2, RuO2_1%, RuO2_2%")
        return
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    
    # x축: 반응 단계 (plot_individual_energetics_diagrams와 동일한 스타일)
    steps = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]
    
    # 색상 설정
    colors = ['black', 'red', 'orange', 'green']
    zorders = [4, 3, 2, 1]
    
    all_energies = []
    
    for i, path_data in enumerate(filtered_data):
        # 에너지 계산
        energy0 = 0.0
        energy1 = path_data['step1']
        energy2 = path_data['step2'] + energy1
        energy3 = path_data['step3'] + energy2
        energy4 = path_data['step4'] + energy3
        energy5 = path_data['step5'] + energy4
        
        # y축: step energies
        step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4, energy5, energy5]
        all_energies.extend(step_energies)
        surface = names[i]
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'], path_data['step5'])
        overpotential = max_energy - 1.23
        # 선 그래프로 그리기 (plot_individual_energetics_diagrams와 동일한 스타일)
        ax.plot(steps, step_energies, '-', color=colors[i % len(colors)], linewidth=2, 
                label=f'{surface} (η = {overpotential:.2f} V)', zorder=zorders[i % len(colors)])
    
    # 그래프 설정 (plot_individual_energetics_diagrams와 동일한 스타일)
    ax.set_xticks([])
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
    
    # y_min = min(all_energies) - 0.2
    # y_max = max(all_energies) + 0.3
    # y_min, y_max = 0.0, 3.0
    # ax.set_ylim(y_min, y_max)
    ax.set_xlim(min(steps), max(steps))
    
    # 최대 에너지와 과전위 정보 추가 (plot_individual_energetics_diagrams와 동일한 스타일)
    info_text = ""
    for path_data in filtered_data:
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'], path_data['step5'])
        overpotential = max_energy - 1.23
        info_text += f'{path_data["surface"]}:\n'
        info_text += f'ΔG1: {path_data["step1"]:.2f} eV\n'
        info_text += f'ΔG2: {path_data["step2"]:.2f} eV\n'
        info_text += f'ΔG3: {path_data["step3"]:.2f} eV\n'
        info_text += f'ΔG4: {path_data["step4"]:.2f} eV\n'
        info_text += f'ΔG5: {path_data["step5"]:.2f} eV\n'
        info_text += f'Overpotential: {overpotential:.2f} V\n\n'
    
    # ax.text(0.02, 0.97, info_text, transform=ax.transAxes, fontsize=10, verticalalignment='top',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.legend(fontsize=12, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    plt.tight_layout()
    # plt.show()
    plt.savefig('all_pathways_combined_5steps_RuO2.png', dpi=300, bbox_inches='tight')
    print("Saved: all_pathways_combined_5steps_RuO2.png")

def plot_all_pathways_combined_ReRuO2(all_paths, oer_data):
    """RuO2, ReRuO2, RuO2_1%, RuO2_2% 케이스를 하나의 그래프에 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting combined pathways")
        return
    
    # 필터링할 surface 목록
    target_surfaces = ['RuO2(d)', 'ReRuO2', 'ReRuO2_1%', 'ReRuO2_2%']
    names = ['RuO$_2$', 'Re(brg)-RuO$_2$', 'Re(brg)-RuO$_2$(+1%)', 'Re(brg)-RuO$_2$(+2%)'] #'Re(brg)-RuO$_2$(+1%)', 'Re(brg)-RuO$_2$(+2%)']

    # 필터링된 데이터만 추출
    filtered_data = []
    
    for path_data in oer_data:
        if path_data['surface'] in target_surfaces:
            filtered_data.append(path_data)
    
    if not filtered_data:
        print("No matching data found for RuO2, ReRuO2, ReRuO2_1%, ReRuO2_2%")
        return
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    
    # x축: 반응 단계 (plot_individual_energetics_diagrams와 동일한 스타일)
    steps = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]
    
    # 색상 설정
    colors = ['black', 'red', 'orange', 'green']
    zorders = [4, 3, 2, 1]
    
    all_energies = []
    
    for i, path_data in enumerate(filtered_data):
        # 에너지 계산
        energy0 = 0.0
        energy1 = path_data['step1']
        energy2 = path_data['step2'] + energy1
        energy3 = path_data['step3'] + energy2
        energy4 = path_data['step4'] + energy3
        energy5 = path_data['step5'] + energy4
        
        # y축: step energies
        step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4, energy5, energy5]
        all_energies.extend(step_energies)
        surface = names[i]
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'], path_data['step5'])
        overpotential = max_energy - 1.23
        # 선 그래프로 그리기 (plot_individual_energetics_diagrams와 동일한 스타일)
        ax.plot(steps, step_energies, '-', color=colors[i % len(colors)], linewidth=2, 
                label=f'{surface} (η = {overpotential:.2f} V)', zorder=zorders[i % len(colors)])
    
    # 그래프 설정 (plot_individual_energetics_diagrams와 동일한 스타일)
    ax.set_xticks([])
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
    
    # y_min = min(all_energies) - 0.2
    # y_max = max(all_energies) + 0.3
    # y_min, y_max = 0.0, 3.0
    # ax.set_ylim(y_min, y_max)
    ax.set_xlim(min(steps), max(steps))
    
    # 최대 에너지와 과전위 정보 추가 (plot_individual_energetics_diagrams와 동일한 스타일)
    info_text = ""
    for path_data in filtered_data:
        max_energy = max(path_data['step1'], path_data['step2'], path_data['step3'], path_data['step4'], path_data['step5'])
        overpotential = max_energy - 1.23
        info_text += f'{path_data["surface"]}:\n'
        info_text += f'ΔG1: {path_data["step1"]:.2f} eV\n'
        info_text += f'ΔG2: {path_data["step2"]:.2f} eV\n'
        info_text += f'ΔG3: {path_data["step3"]:.2f} eV\n'
        info_text += f'ΔG4: {path_data["step4"]:.2f} eV\n'
        info_text += f'ΔG5: {path_data["step5"]:.2f} eV\n'
        info_text += f'Overpotential: {overpotential:.2f} V\n\n'
    
    # ax.text(0.02, 0.97, info_text, transform=ax.transAxes, fontsize=10, verticalalignment='top',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.legend(fontsize=12, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    plt.tight_layout()
    # plt.show()
    plt.savefig('all_pathways_combined_5steps_ReRuO2.png', dpi=300, bbox_inches='tight')
    print("Saved: all_pathways_combined_5steps_ReRuO2.png")

def plot_ruo2_pathways(all_paths, oer_data):
    """RuO2 표면의 Path 1과 Path 2를 하나의 그래프에 그리는 함수"""
    if not all_paths or not oer_data:
        print("No data found for plotting RuO2 pathways")
        return
    
    # Path 1과 Path 2 데이터 찾기
    path1_data = None
    path2_data = None
    
    for i, data in enumerate(oer_data):
        if data['surface'] == 'RuO2(a)' and data['int4'] == 'O_OOH':
            path1_data = data
        elif data['surface'] == 'RuO2(d)' and data['int4'] == 'OH_OO':
            path2_data = data
    
    if not path1_data or not path2_data:
        print("RuO2(a) and RuO2(d) pathways not found in data")
        return
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    
    # x축: 반응 단계 (plot_individual_energetics_diagrams와 동일한 스타일)
    steps = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]

    # 최대 에너지와 과전위 정보 추가
    max1 = max(path1_data['step1'], path1_data['step2'], path1_data['step3'], path1_data['step4'], path1_data['step5'])
    max2 = max(path2_data['step1'], path2_data['step2'], path2_data['step3'], path2_data['step4'], path2_data['step5'])
    over1 = max1 - 1.23
    over2 = max2 - 1.23

    # Path 1: O_V → O_OH → O_O → O_OOH
    energy0 = 0.0
    energy1 = path1_data['step1']
    energy2 = path1_data['step2'] + energy1
    energy3 = path1_data['step3'] + energy2
    energy4 = path1_data['step4'] + energy3
    energy5 = path1_data['step5'] + energy4
    path1_step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4, energy5, energy5]
    ax.plot(steps, path1_step_energies, '-', color='black', linewidth=2, label=f'Single-site mechanism, η = {over1:.2f} V', zorder=2)
    
    # Path 2: O_V → O_OH → O_O → OH_OO
    energy0 = 0.0
    energy1 = path2_data['step1']
    energy2 = path2_data['step2'] + energy1
    energy3 = path2_data['step3'] + energy2
    energy4 = path2_data['step4'] + energy3
    energy5 = path2_data['step5'] + energy4
    path2_step_energies = [energy0, energy0, energy1, energy1, energy2, energy2, energy3, energy3, energy4, energy4, energy5, energy5]
    ax.plot(steps, path2_step_energies, '-', color='red', linewidth=2, label=f'Dual-site mechanism, η = {over2:.2f} V', zorder=1)
    
    # 그래프 설정 (plot_individual_energetics_diagrams와 동일한 스타일)
    ax.set_xticks([])
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
    
    # y축 범위 설정
    all_energies = path1_step_energies + path2_step_energies
    y_min = min(all_energies) - 0.2
    y_max = max(all_energies) + 0.3
    ax.set_ylim(-0.2, 5.4)
    ax.set_xlim(0, 6)
    
    # ax.text(0.02, 0.97, f'Ueff: 2.50 eV\nΔG1: {path1_data["step1"]:.2f} eV\nΔG2: {path1_data["step2"]:.2f} eV\nΔG3a: {path1_data["step3"]:.2f} eV\nΔG3d: {path2_data["step3"]:.2f} eV\nΔG4a: {path1_data["step4"]:.2f} eV\nΔG4d: {path2_data["step4"]:.2f} eV', 
    #         transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # plt.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0), fontsize=10)
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=12)
    plt.tight_layout()
    plt.savefig('RuO2_pathways_comparison.png', dpi=300, bbox_inches='tight')
    # plt.show()
    print("Saved: RuO2_pathways_comparison.png")

if __name__ == "__main__":
    all_paths, oer_data = calculate_oer_energies()
    
    if not all_paths:
        print("OER 에너지 계산에 실패했습니다. 경로와 파일을 확인하세요.")
    else:
        # 표 형태로 데이터 출력
        print_oer_table(oer_data)

        # 모든 pathway를 하나의 그래프에 출력
        print("\nGenerating combined pathways diagram...")
        plot_all_pathways_combined_RuO2(all_paths, oer_data)
        plot_all_pathways_combined_ReRuO2(all_paths, oer_data)
        
        # 개별 energetics diagram 출력
        print("\nGenerating individual energetics diagrams...")
        plot_individual_energetics_diagrams(all_paths, oer_data)

        # RuO2 표면의 Path 1과 Path 2를 하나의 그래프에 그리기
        print("\nGenerating RuO2 pathways comparison...")
        plot_ruo2_pathways(all_paths, oer_data)