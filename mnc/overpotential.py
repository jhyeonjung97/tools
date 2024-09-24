steps = ['H2O->*OH','*OH->*O', '*O->*OOH', '*OOH->O2']
metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']

path = '/pscratch/sd/j/jiuy97/6_MNC/figure'
for metal in metals:
    tsv_file = os.path.join(path, f'{i}_', 'final_with_calculator.json')
    path_pattern = f'/pscratch/sd/j/jiuy97/6_MNC/0_clean/{row_key}/*_{metal}/*_{spin}'
    matching_paths = glob.glob(path_pattern)