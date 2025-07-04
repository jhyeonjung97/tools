from ase.io import read
import os
import pandas as pd

root = '/Users/hailey/Desktop/4_IrFe3'

df = pd.DataFrame(columns=['Ads', 'System', 'Site', 'Ir', 'Fe', 'O'])

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

def get_energy(dir, name):
    atoms_path = os.path.join(root, 'dissolution', dir, 'final.json')
    atoms = read(atoms_path)
    energy = atoms.info['energy']
    vib_path = os.path.join(root, 'dissolution', dir, 'vib.txt')
    vib_correction = read_vib_correction(vib_path)
    results.append({
        'Ads': adsorbate,
        'System': system_name,
        'Site': info['site'],
        **layer_charges,
        **layer_distances,
        **layer_z_avg
    })

dirs = ['1_Ir', '2_IrO_top', '3_IrO_hol', '4_IrFe_Ir', '5_IrFe_Fe', '6_IrFeO_topIr', '7_IrFeO_topFe', '8_IrFeO_holIr', '9_IrFeO_holFe']


print(results)