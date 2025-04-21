import time
import subprocess
import numpy as np
from os import path
from ase.io import read, write
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator
import lightgbm as lgb
from sklearn.model_selection import KFold
import pandas as pd
import argparse

name = 'opt_bulk3_afm'
start_time = time.time()

spin_states_plus_2 = {'Sc': 1, 'Ti': 2, 'V': 3, 'Cr': 4, 'Mn': 5, 'Fe': 4,
                      'Co': 3, 'Ni': 2, 'Cu': 1, 'Zn': 1, 'Ga': 1, 'Ge': 2,
                      'Y': 1, 'Zr': 2, 'Nb': 3, 'Mo': 4, 'Tc': 5, 'Ru': 4,
                      'Rh': 3, 'Pd': 2, 'Ag': 1, 'Cd': 1, 'In': 1, 'Sn': 2,
                      'La': 1, 'Hf': 2, 'Ta': 3, 'W': 4, 'Re': 5, 'Os': 4,
                      'Ir': 3, 'Pt': 2, 'Au': 1, 'Hg': 1, 'Tl': 1, 'Pb': 2
                      }

ldau_luj = {'Ti':{'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr':{'L':2, 'U':3.50,  'J':0.0},
            'Mn':{'L':2, 'U':3.75, 'J':0.0},
            'Fe':{'L':2, 'U':4.30,  'J':0.0},
            'Co':{'L':2, 'U':3.32, 'J':0.0},
            'Ni':{'L':2, 'U':6.45, 'J':0.0},
            'Cu':{'L':2, 'U':9.00,  'J':0.0}
            }

if path.exists('restart.json'):
    atoms = read('restart.json')
elif path.exists('start.traj'):
    atoms = read('start.traj')
    for atom in atoms:
        if atom.symbol in spin_states_plus_2:
            if atom.index % 2 == 1: 
                atom.magmom = spin_states_plus_2[atom.symbol]
            else:
                atom.magmom = -spin_states_plus_2[atom.symbol]
else:
    raise ValueError('Neither restart.json nor start.traj file found')

lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}

def get_kpoints(atoms, l=25, bulk=True):
    cell = atoms.get_cell()
    nkx = int(round(l/np.linalg.norm(cell[0]),0))
    nky = int(round(l/np.linalg.norm(cell[1]),0))
    if bulk == True:
        nkz = int(round(l/np.linalg.norm(cell[2]),0))
    else:
        nkz = 1
    return((nkx, nky, nkz))

kpoints = get_kpoints(atoms, l=25, bulk=True)

atoms.calc = vasp_calculator.Vasp(
                    encut=600,
                    xc='PBE',
                    gga='PE',
                    prec='Normal',
                    #inimix=0,
                    #amix=0.05,
                    #bmix=0.0001,
                    #amix_mag=0.05,
                    #bmix_mag=0.0001,
                    kpts=kpoints,
                    kpar=4,
                    npar=16,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    nelm=200,
                    algo='Normal',
                    isif=3,
                    nsw=200,
                    ibrion=2,
                    ediff=1e-6,
                    ediffg=-0.02,
                    lreal='False',
                    lasph=True, 
                    lvtot=False,
                    laechg=True,
                    # isym=0, 
                    ispin=2,
                    lorbit=11,
                    ldau=True,
                    ldautype=2,
                    ldau_luj=ldau_luj,
                    ldauprint=2,
                    lmaxmix=lmaxmix,
                    setups={'base': 'recommended', 'W': '_sv'},
                    # idipol=3,
                    # dipol=(0, 0, 0.5),
                    # ldipol=True
                    nupdown=0
                    )

energy = atoms.get_potential_energy()
print ('Calculation Complete, storing the run + calculator to traj file')

Trajectory(f'final_{name}.traj','w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("="*40 + "\n")

print(f"Execution time logged in 'time.log'.")

model = lgb.LGBMRegressor(
    n_estimators=200,
    learning_rate=0.05,
    max_depth=4,
    min_data_in_leaf=10,
    min_gain_to_split=0.1,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_alpha=0.1,
    reg_lambda=0.1,
    random_state=42
)

kf = KFold(n_splits=5, shuffle=True, random_state=42)

def select_features_by_correlation(correlation_matrix, target_col, threshold=0.7):
    """
    상관관계가 높은 특성들을 제거하는 함수
    
    Args:
        correlation_matrix: 상관관계 행렬
        target_col: 타겟 변수 이름
        threshold: 상관관계 임계값 (0.7, 0.8, 0.9 등)
    
    Returns:
        selected_features: 선택된 특성들의 리스트
    """
    # 타겟 변수와의 상관관계로 정렬
    target_corr = correlation_matrix[target_col].abs().sort_values(ascending=False)
    
    # 선택된 특성과 제거된 특성 추적
    selected_features = []
    dropped_features = []
    
    # 모든 특성에 대해 반복
    for feature in target_corr.index:
        if feature == target_col:
            continue
            
        # 이미 선택된 특성들과의 상관관계 확인
        high_corr = False
        for selected in selected_features:
            if abs(correlation_matrix.loc[feature, selected]) > threshold:
                high_corr = True
                # 타겟과의 상관관계가 더 높은 특성 선택
                if target_corr[feature] > target_corr[selected]:
                    selected_features.remove(selected)
                    dropped_features.append(selected)
                    selected_features.append(feature)
                else:
                    dropped_features.append(feature)
                break
                
        if not high_corr:
            selected_features.append(feature)
    
    print(f"\n상관관계 임계값: {threshold}")
    print("\n선택된 특성:")
    for feat in selected_features:
        print(f"- {feat}")
    
    print("\n제거된 특성:")
    for feat in dropped_features:
        print(f"- {feat}")
    
    return selected_features

# 상관관계 임계값을 파라미터로 받기
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--corr_threshold', type=float, default=0.7,
                   help='상관관계 임계값 (0.7, 0.8, 0.9 등)')

# 상관관계 기반 특성 선택
selected_features = select_features_by_correlation(
    correlation_matrix=df[args.X].corr(),
    target_col=args.Y,
    threshold=args.corr_threshold
)

# 선택된 특성으로 X 데이터 업데이트
X = X[selected_features]