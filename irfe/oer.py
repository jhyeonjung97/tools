import pandas as pd
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import os

# 기본 경로 설정
base_path = '/Users/hailey/Desktop/4_IrFe3'

# 반응 경로 정의
mechanism1 = ['1_V_V', '3_V_OH', '2_V_O', '5_O_OH']  # 첫 번째 메커니즘
mechanism2 = ['2_V_O', '5_O_OH', '4_O_O', '6_O_OOH']    # 두 번째 메커니즘
mechanism3 = ['1_V_V', '3_V_OH', '7_OH_OH', '5_O_OH']  # 세 번째 메커니즘 (LOM')

# 표면 정의
surfaces = ['1_Ir', '2_Ir', '3_Fe', '4_Fe', '5_Fe', '6_Fe']
surface_labels = ['Ir1', 'Ir2', 'IrFe1', 'IrFe2', 'IrFe3', 'IrFe4']

# 결과를 저장할 데이터프레임 생성
results = pd.DataFrame(index=surfaces, columns=['AEM1', 'AEM2', 'LOMa1', 'LOMa2', 'LOMb1', 'LOMb2'])

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

# 각 중간체의 에너지와 dG를 저장할 리스트
energy_data = []

# 각 표면과 메커니즘에 대해 과전압 계산
for surface in surfaces:
    # R1과 R2에 대해 계산
    for r_folder in ['4_R_top', '5_R_hol']:
        # 첫 번째 메커니즘
        energies1 = []
        for step in mechanism1:
            path = f'{base_path}/{r_folder}/{surface}/{step}'
            energy = get_energy(path)
            energies1.append(energy)
        print(energies1)
        # 두 번째 메커니즘
        energies2 = []
        for step in mechanism2:
            path = f'{base_path}/{r_folder}/{surface}/{step}'
            energy = get_energy(path)
            energies2.append(energy)
        print(energies2)
        # 세 번째 메커니즘
        energies3 = []
        for step in mechanism3:
            path = f'{base_path}/{r_folder}/{surface}/{step}'
            energy = get_energy(path)
            energies3.append(energy)
        print(energies3)
        # 과전압 계산
        if None not in energies1:
            overpotential1, dG1_1, dG2_1, dG3_1, dG4_1 = calculate_overpotential(energies1)
            # 첫 번째 메커니즘 데이터 저장
            energy_data.append({
                'Surface': surface,
                'rxn': r_folder,
                'Mechanism': 'LOMa',
                'Step1': mechanism1[0],
                'Step2': mechanism1[1],
                'Step3': mechanism1[2],
                'Step4': mechanism1[3],
                'Energy1': energies1[0],
                'Energy2': energies1[1],
                'Energy3': energies1[2],
                'Energy4': energies1[3],
                'dG1': dG1_1,
                'dG2': dG2_1,
                'dG3': dG3_1,
                'dG4': dG4_1,
                'Overpotential': overpotential1
            })
        
        if None not in energies2:
            overpotential2, dG1_2, dG2_2, dG3_2, dG4_2 = calculate_overpotential(energies2)
            # 두 번째 메커니즘 데이터 저장
            energy_data.append({
                'Surface': surface,
                'rxn': r_folder,
                'Mechanism': 'AEM',
                'Step1': mechanism2[0],
                'Step2': mechanism2[1],
                'Step3': mechanism2[2],
                'Step4': mechanism2[3],
                'Energy1': energies2[0],
                'Energy2': energies2[1],
                'Energy3': energies2[2],
                'Energy4': energies2[3],
                'dG1': dG1_2,
                'dG2': dG2_2,
                'dG3': dG3_2,
                'dG4': dG4_2,
                'Overpotential': overpotential2
            })
            
        if None not in energies3:
            overpotential3, dG1_3, dG2_3, dG3_3, dG4_3 = calculate_overpotential(energies3)
            # 세 번째 메커니즘 데이터 저장
            energy_data.append({
                'Surface': surface,
                'rxn': r_folder,
                'Mechanism': 'LOMb',
                'Step1': mechanism3[0],
                'Step2': mechanism3[1],
                'Step3': mechanism3[2],
                'Step4': mechanism3[3],
                'Energy1': energies3[0],
                'Energy2': energies3[1],
                'Energy3': energies3[2],
                'Energy4': energies3[3],
                'dG1': dG1_3,
                'dG2': dG2_3,
                'dG3': dG3_3,
                'dG4': dG4_3,
                'Overpotential': overpotential3
            })
        
        # 결과 저장
        if r_folder == '4_R_top':
            results.loc[surface, 'LOMa1'] = overpotential1 if None not in energies1 else None
            results.loc[surface, 'AEM1'] = overpotential2 if None not in energies2 else None
            results.loc[surface, 'LOMb1'] = overpotential3 if None not in energies3 else None
        elif r_folder == '5_R_hol':
            results.loc[surface, 'LOMa2'] = overpotential1 if None not in energies1 else None
            results.loc[surface, 'AEM2'] = overpotential2 if None not in energies2 else None
            results.loc[surface, 'LOMb2'] = overpotential3 if None not in energies3 else None

# 에너지 데이터를 데이터프레임으로 변환
energy_df = pd.DataFrame(energy_data)

# float 값을 소수점 두 자리까지만 표시
float_columns = ['Energy1', 'Energy2', 'Energy3', 'Energy4', 'dG1', 'dG2', 'dG3', 'dG4', 'Overpotential']
energy_df[float_columns] = energy_df[float_columns].round(2)

# TSV 파일로 저장
energy_df.to_csv('oer_energies.tsv', sep='\t', index=False)

# 결과 출력
print("\nOER 과전압 (V):")
print(results)

# 결과 시각화
plt.figure(figsize=(4, 3))
ax = results.plot(kind='bar', ax=plt.gca(), 
                #  color=['black', 'red', 'blue', 'green', 'yellow', 'purple'],
                 color=['black', 'white'], edgecolor='black', hatch='////'
                 )
plt.ylabel('OER Overpotential (V)')
plt.xticks(range(len(surface_labels)), surface_labels, rotation=45)
plt.legend(bbox_to_anchor=(0.5, 1.1), loc='center', ncol=4)
plt.tight_layout()
plt.savefig('OER_overpotential.png', bbox_inches='tight')
plt.close()

# 각 표면과 메커니즘별로 꺾은선 그래프 그리기
for surface in surfaces:
    for mech, label in zip(['LOMa', 'AEM', 'LOMb'], ['LOMa', 'AEM', 'LOMb']):
        for rxn, rxn_label in zip(['4_R_top', '5_R_hol'], ['1', '2']):
            df = energy_df[(energy_df['Surface'] == surface) & (energy_df['Mechanism'] == mech) & (energy_df['rxn'] == rxn)]
            if df.empty:
                continue
            dgs = df[['dG1', 'dG2', 'dG3', 'dG4']].values.flatten()
            yy = [0, 0, dgs[0], dgs[0], dgs[0]+dgs[1], dgs[0]+dgs[1], dgs[0]+dgs[1]+dgs[2], dgs[0]+dgs[1]+dgs[2], dgs[0]+dgs[1]+dgs[2]+dgs[3], dgs[0]+dgs[1]+dgs[2]+dgs[3]]
            plt.figure(figsize=(4,3))
            plt.plot(range(10), yy, color='black')
            plt.xlabel('Reaction coordinate')
            plt.ylabel('Relative energy (ΔG, eV)')
            plt.xlim(0, 9)
            plt.ylim(0, 5)
            plt.tight_layout()
            plt.xticks([])  # x축 눈금 제거
            plt.savefig(f'OER_{surface_labels[surfaces.index(surface)]}_{label}{rxn_label}_profile.png', bbox_inches='tight')
            plt.close()