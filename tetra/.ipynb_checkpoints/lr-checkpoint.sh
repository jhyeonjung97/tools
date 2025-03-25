    parser = argparse.ArgumentParser(description='Linear regression using bulk_data.csv and mendeleev_data.csv')
    parser.add_argument('--Y', default='form', help='Target column from bulk_data.csv (default: form)')
    parser.add_argument('--X', nargs='+', default=[
        'chg', 'mag', 'volume', 'l_bond', 'n_bond',
        'grosspop', 'madelung', '-ICOHPm', '-ICOHPn', 'ICOBIm', 'ICOBIn', '-ICOOPm', '-ICOOPn', 
        'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 
        'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 
        'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform',
    ], help='List of feature columns from bulk_data.csv and/or mendeleev_data.csv')
    parser.add_argument('--row', nargs='+', type=str, default=None, help='Filter by row: 3d, 4d, or 5d')
    parser.add_argument('--coord', nargs='+', type=str, default=None, help='Filter by coordination, e.g., ZB, RS')
    parser.add_argument('--output', type=str, default='result', help='Output filename suffix')
    args = parser.parse_args()
    
python ~/bin/tools/tetra/lr.py --Y 'form' --X ['chg', 'mag', 'volume', 'l_bond', 'n_bond', 'grosspop', 'madelung', '-ICOHPm', '-ICOHPn', 'ICOBIm', 'ICOBIn', '-ICOOPm', '-ICOOPn', 'pauling', 'ion1', 'ion2', 'ion12', 'ion3', 'Natom', 'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', 'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform'] --row ['3d', '4d', '5d']