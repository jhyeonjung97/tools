import os
import socket

hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/8_V_slab'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/8_V_slab'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/8_V_slab'
else:
    raise ValueError(f"Unknown hostname: {hostname} or user_name: {user_name}. Please set the root path manually.")
save_path = os.path.join(root, 'figures')

metals = {
    '3d': ['Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge'],
    '4d': ['Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'],
    '5d': ['Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb'],
}

surface_names = ['corundum', 'rutile', 'C12c', 'Pm3m']

for c, coord in enumerate(['+3', '+4', '+5', '+6']):
    surface_name = surface_names[c]
    for row in ['3d', '4d', '5d']:
        for m, metal in enumerate(metals[row]):
            m_str = f"{m:02d}"
            if coord == '+3':
                from_path = os.path.join(root, 'cathub', f'{metal}2O3-corundum', '012')
                to_path = os.path.join(root, 'comer', '1_Octahedral_+3_012', row, f'{m_str}_{metal}')
                os.system(f'cp {from_path}/*.json {from_path}/*/*.json {to_path}')
            elif coord == '+4':
                from_path = os.path.join(root, 'cathub', f'{metal}O2-rutile_rutile', '100')
                to_path = os.path.join(root, 'comer', '2_Octahedral_+4_100', row, f'{m_str}_{metal}')
                os.system(f'cp {from_path}/*.json {from_path}/*/*.json {to_path}')

                from_path = os.path.join(root, 'cathub', f'{metal}O2-rutile_rutile', '110')
                to_path = os.path.join(root, 'comer', '3_Octahedral_+4_110', row, f'{m_str}_{metal}')
                os.system(f'cp {from_path}/*.json {from_path}/*/*.json {to_path}')
            elif coord == '+5':
                from_path = os.path.join(root, 'cathub', f'{metal}2O5-C12c', '100-lc')
                to_path = os.path.join(root, 'comer', '4_Octahedral_+5_100', row, f'{m_str}_{metal}')
                os.system(f'cp {from_path}/*.json {from_path}/*/*.json {to_path}')
            elif coord == '+6':
                from_path = os.path.join(root, 'cathub', f'{metal}O3-Pm3m', '001')
                to_path = os.path.join(root, 'comer', '5_Octahedral_+6_001', row, f'{m_str}_{metal}')
                os.system(f'cp {from_path}/*.json {from_path}/*/*.json {to_path}')

# find . -name 'OH.json' -execdir mv OH.json HO.json \;
