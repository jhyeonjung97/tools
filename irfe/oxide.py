from ase.io import read
import os
import pandas as pd

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

root = '/Users/hailey/Desktop/4_IrFe3/'
surfaces = ['Ir', 'IrFe']
df = pd.DataFrame(columns=surfaces, index=range(1, 9))

for s in range(len(surfaces)):
    clean_path = os.path.join(root, 'slab', f'{s*6}_{surfaces[s]}')
    atoms_path = os.path.join(clean_path, 'final.json')
    vib_path = os.path.join(clean_path, 'vib.txt')
    energy = read(atoms_path).get_potential_energy()
    vib_correction = read_vib_correction(vib_path)
    clean_energy = energy + vib_correction
    
    ads_path = os.path.join(root, '3_O', f'{s*6}_{surfaces[s]}', '1_layer_top')
    # if s == 0:
    #     ads_path = os.path.join(root, '3_O', f'{s*6}_{surfaces[s]}', '2_layer_hol')
    atoms_path = os.path.join(ads_path, 'final.json')
    vib_path = os.path.join(ads_path, 'vib-z', 'vib.txt')
    energy = read(atoms_path).get_potential_energy()
    vib_correction = read_vib_correction(vib_path)
    ads_energy = energy + vib_correction

    for i in range(1, 9):
        oxide_path = os.path.join(root, 'oxide', f'{s+7}_{surfaces[s]}Ox', f'{i}_')
        atoms_path = os.path.join(oxide_path, 'final.json')
        vib_path = os.path.join(oxide_path, 'vib.txt')
        energy = read(atoms_path).get_potential_energy()
        vib_correction = read_vib_correction(vib_path)
        oxide_energy = energy + vib_correction

        potential = ((oxide_energy - clean_energy)/16 - (G_H2O - G_H2)) / 2
        energy_diff = oxide_energy - clean_energy
        # df.loc[i, surfaces[s]] = potential
        df.loc[i, surfaces[s]] = oxide_energy

for i in range(1, 9):
    df.loc[i, 'diff'] = df.loc[i, 'Ir'] - df.loc[i, 'IrFe']

        potential = ((oxide_energy - ads_energy)/8 - (G_H2O - G_H2)) / 2
        df.loc[i, f'{surfaces[s]}O'] = potential

print(df)
df.to_csv('/Users/hailey/Desktop/4_IrFe3/figures/oxide.csv')