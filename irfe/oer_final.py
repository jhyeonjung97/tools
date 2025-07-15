import pandas as pd
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import os

# 기본 경로 설정
base_path = '/Users/hailey/Desktop/4_IrFe3'

# 반응 경로 정의 (LOMa 메커니즘만 사용)
mechanism = ['1_V_V', '3_V_OH', '2_V_O', '5_O_OH']

# 표면 정의 (두 개만 사용)
surfaces = ['1_Ir_top', '3_IrFe_top1']
surface_labels = ['Ir', 'IrFe']

# 기체 분자의 에너지와 보정값
E_H2 = -6.77108058  # H2의 에너지 (eV)
E_H2O = -14.22266334  # H2O의 에너지 (eV)

# ZPE, 열용량, 엔트로피 보정
ZPE_H2O = 0.558
CV_H2O = 0.103
TS_H2O = 0.675

ZPE_H2 = 0.268
CV_H2 = 0.0905
TS_H2 = 0.408

# 기체 분자의 깁스 자유에너지
G_H2O = E_H2O + ZPE_H2O - TS_H2O + CV_H2O
G_H2 = E_H2 + ZPE_H2 - TS_H2 + CV_H2

def count_atoms(atoms):
    """원자 개수를 세는 함수"""
    symbols = atoms.get_chemical_symbols()
    n_H = symbols.count('H')
    n_O = symbols.count('O')
    return n_H, n_O

def get_energy(path):
    """에너지 계산 함수"""
    # 필요한 파일 체크
    if not os.path.exists(f'{path}/DONE'):
        return None
    if not os.path.exists(f'{path}/vib.txt'):
        return None
    if os.path.exists(f'{path}/unmatched'):
        return None
        
    try:
        atoms = read(f'{path}/final.json')
        energy = atoms.get_potential_energy()
        
        # 진동 계산 결과가 있는 경우 ZPE와 TS 보정
        try:
            with open(f'{path}/vib.txt', 'r') as f:
                content = f.read()
                zpe = 0
                ts = 0
                for line in content.split('\n'):
                    if 'Zero-point energy E_ZPE' in line:
                        zpe = float(line.split()[-2])
                    if 'Entropy contribution T*S' in line:
                        ts = float(line.split()[-2])
                energy += zpe - ts
        except:
            pass
            
        # 화학양론 보정
        n_H, n_O = count_atoms(atoms)
        
        # n_H * (G_H2 / 2) + n_O * (G_H2O - G_H2)로 보정
        energy -= n_H * (G_H2 / 2) + n_O * (G_H2O - G_H2)
            
        return energy
    except:
        return None

def calculate_overpotential(energies):
    """과전압 계산 함수"""
    if None in energies:
        return None
    
    # 각 단계의 깁스 자유에너지 변화 계산
    dG1 = energies[1] - energies[0]  # *O+* -> *O+*OH or *+* -> *+*OH
    dG2 = energies[2] - energies[1]  # *O+*OH -> *O+*O or *+*OH -> *+*O
    dG3 = energies[3] - energies[2]  # *O+*O -> *O+*OOH or *+*O -> *OH+*O
    dG4 = 4.92 - (dG1 + dG2 + dG3)  # 마지막 단계
    
    # 최대 깁스 자유에너지 변화 찾기
    max_dG = max(dG1, dG2, dG3, dG4)
    
    # 과전압 계산 (1.23V가 이론값)
    overpotential = max_dG - 1.23
    
    return overpotential, dG1, dG2, dG3, dG4

# 에너지 데이터를 저장할 리스트
energy_data = []

# 각 표면에 대해 LOMa2 메커니즘 계산
for surface in surfaces:
    # R2에 대해서만 계산 (LOMa2)
    r_folder = '5_R2'
    
    energies = []
    for step in mechanism:
        path = f'{base_path}/{r_folder}/{surface}/{step}'
        energy = get_energy(path)
        energies.append(energy)
    
    # 과전압 계산
    if None not in energies:
        overpotential, dG1, dG2, dG3, dG4 = calculate_overpotential(energies)
        # 데이터 저장
        energy_data.append({
            'Surface': surface,
            'Surface_Label': surface_labels[surfaces.index(surface)],
            'Mechanism': 'LOMa2',
            'Step1': mechanism[0],
            'Step2': mechanism[1],
            'Step3': mechanism[2],
            'Step4': mechanism[3],
            'Energy1': energies[0],
            'Energy2': energies[1],
            'Energy3': energies[2],
            'Energy4': energies[3],
            'dG1': dG1,
            'dG2': dG2,
            'dG3': dG3,
            'dG4': dG4,
            'Overpotential': overpotential
        })

# 에너지 데이터를 데이터프레임으로 변환
energy_df = pd.DataFrame(energy_data)

# float 값을 소수점 두 자리까지만 표시
float_columns = ['Energy1', 'Energy2', 'Energy3', 'Energy4', 'dG1', 'dG2', 'dG3', 'dG4', 'Overpotential']
energy_df[float_columns] = energy_df[float_columns].round(2)

# TSV 파일로 저장
energy_df.to_csv('oer_final_energies.tsv', sep='\t', index=False)

# 결과 출력
print("\nOER 과전압 (V) - LOMa2 메커니즘:")
for _, row in energy_df.iterrows():
    print(f"{row['Surface_Label']}: {row['Overpotential']:.2f} V")

# 두 에너지 프로파일을 하나의 이미지에 그리기
plt.figure(figsize=(4, 3))

for i, (_, row) in enumerate(energy_df.iterrows()):
    dgs = [row['dG1'], row['dG2'], row['dG3'], row['dG4']]
    
    # 원본 oer.py와 동일한 방식으로 에너지 프로파일 계산
    yy = [0, 0, dgs[0], dgs[0], dgs[0]+dgs[1], dgs[0]+dgs[1], dgs[0]+dgs[1]+dgs[2], dgs[0]+dgs[1]+dgs[2], dgs[0]+dgs[1]+dgs[2]+dgs[3], dgs[0]+dgs[1]+dgs[2]+dgs[3]]
    
    # 반응 좌표 (0부터 9까지)
    reaction_coord = list(range(10))
    
    if i == 0:
        plt.plot(reaction_coord, yy, 
                color='black', marker=None, dashes=[3, 1],
                label=f"{row['Surface_Label']} (η = {row['Overpotential']:.2f} V)")
    else:
        plt.plot(reaction_coord, yy, 
                color='black', marker=None,
                label=f"{row['Surface_Label']} (η = {row['Overpotential']:.2f} V)")

plt.xlabel('Reaction Coordinate')
plt.ylabel('Relative Energy (ΔG, eV)')
plt.legend(loc='upper left', handletextpad=0.5, labelspacing=0.2)
plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5], ['*+H$_2$O', '*OH', '*O', '*OOH', r'*+O$_2$'])
plt.xlim(0, 9)
plt.ylim(0, 6)
plt.tight_layout()

# figures 폴더가 없으면 생성
os.makedirs(f'{base_path}/figures', exist_ok=True)

plt.savefig(f'{base_path}/figures/OER_final_profiles.png', 
            bbox_inches='tight', transparent=True, dpi=300)
plt.show()
plt.close()

print(f"\n에너지 프로파일이 {base_path}/figures/OER_final_profiles.png에 저장되었습니다.") 