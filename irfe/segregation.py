import pandas as pd
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt


adsorbates = ['None', 'H_1ML', 'O_1ML']
tm_ads_dirs = ['0_slab', '2_H_hol', '4_O_hol'] # '3_O_top',
surface_dirs = ['5_IrMn', '6_IrFe', '7_IrCo', '8_IrNi']
surfaces = ['IrMn', 'IrFe', 'IrCo', 'IrNi']

df = pd.DataFrame()

for i, tm_ads_dir in enumerate(tm_ads_dirs):
    for j, surface_dir in enumerate(surface_dirs):
        adsorbate = adsorbates[i]
        surface = surfaces[j]

        tm_path = f'/Users/hailey/Desktop/4_IrFe3/sgr/{tm_ads_dir}/{surface_dir}'
        tm_atom_path = f'{tm_path}/final.json'
        tm_atoms = read(tm_atom_path)
        tm_energy = tm_atoms.get_potential_energy()

        if tm_ads_dir == '0_slab':
            ir_path = f'/Users/hailey/Desktop/4_IrFe3/slab/{surface_dir}'
        elif tm_ads_dir == '2_H_hol':
            ir_path = f'/Users/hailey/Desktop/4_IrFe3/1_H/{surface_dir}/1_layer_top'
        # elif tm_ads_dir == '3_O_top':
        #     ir_atom_path = f'/Users/hailey/Desktop/4_IrFe3/3_O/{surface_dir}/1_layer_top'
        elif tm_ads_dir == '4_O_hol':
            ir_path = f'/Users/hailey/Desktop/4_IrFe3/3_O/{surface_dir}/1_layer_top'
        ir_atom_path = f'{ir_path}/final.json'
        ir_atoms = read(ir_atom_path)
        ir_energy = ir_atoms.get_potential_energy()

        if tm_ads_dir != '0_slab':

            tm_vib_path = f'{tm_path}/vib.txt'
            with open(tm_vib_path, 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if 'Zero-point energy E_ZPE' in line:
                        zpe = float(line.split()[-2])
                    if 'Entropy contribution T*S' in line:
                        ts = float(line.split()[-2])
                print(zpe, ts)
                tm_energy += zpe - ts
                
            ir_vib_path = f'{ir_path}/vib.txt'
            with open(ir_vib_path, 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if 'Zero-point energy E_ZPE' in line:
                        zpe = float(line.split()[-2])
                    if 'Entropy contribution T*S' in line:
                        ts = float(line.split()[-2])
                ir_energy += zpe - ts

        cell = ir_atoms.get_cell()
        surface_area = np.linalg.norm(np.cross(cell[0], cell[1]))
        
        df.loc[surface, adsorbate] = (ir_energy - tm_energy) / surface_area

print(df)

# Create the plot
plt.figure(figsize=(4, 3))
df.plot(marker='o', ax=plt.gca())
plt.axhline(y=0, color='black', linewidth=1.0)
plt.ylim(-0.1, 0.4)
plt.ylabel('ΔG / Surface Area (eV/Å²)')
plt.legend()
plt.tight_layout()
plt.savefig('segregation_plot.png')
plt.show()
plt.close()
