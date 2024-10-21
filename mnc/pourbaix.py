rows = {
    '3d': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu'],
    '4d': ['Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd'],
    '5d': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']
}

for row in rows:
    for m, metal in enumerate(rows[row]):
        i = m + 2
        for spin in {'LS', 'IS', 'HS'}:
            for dz in {'1', '5'}:
                dir_ = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row}/{i}_{metal}/*_{spin}/{dz}_'
                dir_o = f'/pscratch/sd/j/jiuy97/6_MNC/1_O/{row}/{i}_{metal}/*_{spin}/{dz}_'
                dir_oh = f'/pscratch/sd/j/jiuy97/6_MNC/1_O/{row}/{i}_{metal}/*_{spin}/{dz}_'
                atoms = 
        
        