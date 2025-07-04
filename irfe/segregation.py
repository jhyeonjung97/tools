from ase.io import read
import os
import pandas as pd

root = '/Users/hailey/Desktop/4_IrFe3'

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

df = pd.DataFrame(columns=['metal', 'o-covered', 'oxide-l1', 'oxide-l2'], index=['Ir', 'Ir-Fe'])

def read_vib_correction(vib_file):
    try:
        with open(vib_file, 'r') as f:
            content = f.read()
            for line in content.split('\n'):
                if 'Zero-point energy E_ZPE' in line:
                    zpe = float(line.split()[-2])
                if 'Entropy contribution T*S' in line:
                    ts = float(line.split()[-2])
            return zpe - ts
    except:
        return 0.0

def get_energy(dir, row, col):
    atoms_path = os.path.join(root, 'oxide/9_two_layers', dir, 'final.json')
    atoms = read(atoms_path)
    energy = atoms.get_potential_energy()
    vib_path = os.path.join(root, 'oxide/9_two_layers', dir, 'vib.txt')
    if dir == '4_IrIrOIrO':
        vib_path = os.path.join(root, 'oxide/9_two_layers', '8_IrFeIrOIrO', 'vib.txt')
    vib_correction = read_vib_correction(vib_path)
    df.loc[row, col] = energy + vib_correction

dirs = ['1_IrIrIr', '2_IrIrIr', '3_IrIrIrO', '4_IrIrOIrO', '5_IrFeIrIr', '6_IrFeIrIr', '7_IrFeIrIrO', '8_IrFeIrOIrO']

get_energy(dirs[0], 'Ir', 'metal')
get_energy(dirs[1], 'Ir', 'o-covered')
get_energy(dirs[2], 'Ir', 'oxide-l1')
get_energy(dirs[3], 'Ir', 'oxide-l2')
get_energy(dirs[4], 'Ir-Fe', 'metal')
get_energy(dirs[5], 'Ir-Fe', 'o-covered')
get_energy(dirs[6], 'Ir-Fe', 'oxide-l1')
get_energy(dirs[7], 'Ir-Fe', 'oxide-l2')

df['dg(metal-*O)'] = ((df['o-covered'] - df['metal']) / 8 - (G_H2O - G_H2)) / 2
df['dg(*O-oxide)'] = ((df['oxide-l1'] - df['o-covered']) / 8 - (G_H2O - G_H2)) / 2
df['dg(l1-l2)'] = ((df['oxide-l2'] - df['oxide-l1']) / 16 - (G_H2O - G_H2)) / 2
df['dg(l0-l1)'] = ((df['oxide-l1'] - df['metal']) / 16 - (G_H2O - G_H2)) / 2
df['dg24'] = ((df['oxide-l2'] - df['o-covered']) / 24 - (G_H2O - G_H2)) / 2
df['dg14'] = ((df['oxide-l2'] - df['metal']) / 32 - (G_H2O - G_H2)) / 2

print(df.to_string(float_format=lambda x: f"{x:.2f}"))