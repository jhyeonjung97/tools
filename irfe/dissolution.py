from ase.io import read
import os
import pandas as pd

results = []

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

def get_energy(dir):
    atoms_path = os.path.join(dir, 'final.json')
    atoms = read(atoms_path)
    energy = atoms.get_potential_energy()
    vib_path = os.path.join(dir, 'vib-z.txt')
    vib_correction = read_vib_correction(vib_path)
    return energy, vib_correction

dirs = {'1_Ir': 'slab/0_Ir',
        '4_IrFe_Ir': 'slab/2_Fe',
        '5_IrFe_Fe': 'slab/2_Fe',
        '2_IrO_top': '3_O/0_Ir/1_layer_top',
        '3_IrO_hol': '3_O/0_Ir/2_layer_hol',
        '6_IrFeO_topIr': '3_O/2_Fe/1_layer_top',
        '7_IrFeO_topFe': '3_O/2_Fe/1_layer_top',
        '8_IrFeO_holIr': '3_O/2_Fe/2_layer_hol',
        '9_IrFeO_holFe': '3_O/2_Fe/2_layer_hol'}

for dir in dirs.keys():
    system_name = dir.split('_')[1]

    if 'top' in dir:
        site = 'top'
    elif 'hol' in dir:
        site = 'hol'
    else:
        site = '-'

    before_path = os.path.join('/Users/hailey/Desktop/4_IrFe3', dirs[dir])
    after_path = os.path.join('/Users/hailey/Desktop/4_IrFe3/dissolution', dir)
    energy_before_dissolution = get_energy(before_path)
    energy_after_dissolution = get_energy(after_path)

    if dir.endswith('Fe'):
        metal_name = 'Fe'
        electron = 2
        metal_energy = -32.32542329/4
        standard_potential = -0.44
    else:
        metal_name = 'Ir'
        electron = 3
        metal_energy = -35.36938889/4
        standard_potential = 1.156
    ion_energy = metal_energy + electron * standard_potential

    energy_before_dissolution, vib_before_dissolution = energy_before_dissolution
    energy_after_dissolution, vib_after_dissolution = energy_after_dissolution
    gibbs_before = energy_before_dissolution + vib_before_dissolution
    gibbs_after = energy_after_dissolution + vib_after_dissolution + ion_energy
    potential = (gibbs_after - gibbs_before) / electron



    results.append({
        'System': system_name,
        'Site': site,
        'Metal': metal_name,
        'Ion_energy': ion_energy,
        'E_before': energy_before_dissolution,
        'E_after': energy_after_dissolution,
        'V_before': vib_before_dissolution,
        'V_after': vib_after_dissolution,
        'Potential': potential,
    })

df = pd.DataFrame(results)
print(df)
# df.to_csv('dissolution.csv', index=False)