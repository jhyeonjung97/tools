import argparse
import glob
import json
import os
import re
from collections import Counter
from math import log10

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.io import read, write
from mendeleev import element
from pymatgen.core.ion import Ion

SUBSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉',}
SUPERSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',}
DIVERGING_COLORMAPS = ['RdBu', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',]

# Thermodynamic constants (set by init_thermo_constants)
const = water_formation_energy = gh = go = goh = dgh = dgo = dgoh = None

def init_thermo_constants():
    """Initialize module-level thermodynamic constants."""
    global const, water_formation_energy, gh, go, goh, dgh, dgo, dgoh

    calmol = 23.061
    kb = 8.617e-5  # eV/K
    T = 298.15  # K
    const = kb * T * np.log(10)  # 0.0592 eV
    water_formation_energy = 56.690 / calmol

    h2 = -6.77149190
    h2o = -14.23091949
    gh2o = h2o + 0.558 - 0.675 + 0.103
    gh2 = h2 + 0.268 - 0.408 + 0.0905

    gh = gh2 / 2
    go = gh2o - gh2
    goh = gh2o - gh2 / 2

    dgo = 0.064 + 0.034 - 0.060
    dgoh = 0.376 + 0.042 - 0.066
    dgh = dgoh - dgo

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')
    parser.add_argument('--json-dir', type=str, default='.', help='Folder with json files')
    parser.add_argument('--csv-dir', type=str, default='.', help='Folder with label.csv')
    parser.add_argument('--suffix', type=str, default='', help='Output filename suffix')
    parser.add_argument('--hybrid', action='store_true', help='Hybrid mode')
    parser.add_argument('--gc', action='store_true', help='Grand Canonical DFT (A, B, C columns)')
    parser.add_argument('--pH', type=int, default=0, help='pH for 1D plot (default: 0)')
    parser.add_argument('--concentration', type=float, default=1e-6, help='Ion concentration (M)')
    parser.add_argument('--pressure', type=float, default=1e-6, help='Gas pressure (atm)')
    parser.add_argument('--gibbs', action='store_true', help='Apply G_corr from label.csv')
    parser.add_argument('--oh-min', type=float, default=0.0,
                        help='Min O-H distance (Å) for bond counting (VESTA default: 0)')
    parser.add_argument('--oh-max', type=float, default=1.1,
                        help='Max O-H distance (Å) for bond counting (VESTA default: 1.1)')
    parser.add_argument('--tick', type=float, default=0.01, help='Grid tick size')
    parser.add_argument('--pHmin', type=float, default=0)
    parser.add_argument('--pHmax', type=float, default=14)
    parser.add_argument('--Umin', type=float, default=-1)
    parser.add_argument('--Umax', type=float, default=3)
    parser.add_argument('--Gmin', type=float, help='Y-axis min for 1D plot')
    parser.add_argument('--Gmax', type=float, help='Y-axis max for 1D plot')
    parser.add_argument('--figx', type=float, default=4)
    parser.add_argument('--figy', type=float, default=3)
    parser.add_argument('--HER', action='store_true')
    parser.add_argument('--OER', action='store_true')
    parser.add_argument('--legend-in', action='store_true')
    parser.add_argument('--legend-out', action='store_true')
    parser.add_argument('--legend-up', action='store_true')
    parser.add_argument('--cmap', type=str, default='Greys', help='Colormap for bulk/combination')
    parser.add_argument('--cmin', type=float, default=0.1)
    parser.add_argument('--cmax', type=float, default=0.7)
    parser.add_argument('--cgap', type=float, default=0.0)
    parser.add_argument('--cmap-2d', type=str, default='RdBu', help='Colormap for original surfaces')
    parser.add_argument('--cmin-2d', type=float, default=0.0)
    parser.add_argument('--cmax-2d', type=float, default=1.0)
    parser.add_argument('--cgap-2d', type=float, default=0.2)
    parser.add_argument('--cmap-1d', type=str, default='Spectral')
    parser.add_argument('--cmin-1d', type=float, default=0.0)
    parser.add_argument('--cmax-1d', type=float, default=1.0)
    parser.add_argument('--cgap-1d', type=float, default=0.0)
    parser.add_argument('--no-bulk', action='store_true')
    parser.add_argument('--show-fig', action='store_true')
    parser.add_argument('--show-thermo', action='store_true')
    parser.add_argument('--show-ref', action='store_true')
    parser.add_argument('--show-element', action='store_true')
    parser.add_argument('--show-count', action='store_true')
    parser.add_argument('--show-label', action='store_true')
    parser.add_argument('--show-oh', action='store_true', help='Show auto-detected O-H bond counts')
    parser.add_argument('--show-min-coord', action='store_true')
    parser.add_argument('--png', action='store_true')
    parser.add_argument('--png-rotation', type=str, default='-90x, -90y, 0z')
    parser.add_argument('--label-csv', type=str)
    parser.add_argument('--thermo-data', type=str)
    parser.add_argument('--ref-energies', type=str)
    return parser.parse_args()

def build_output_names(args):
    """Build output PNG basename and suffix from CLI args."""
    png_name = 'pourbaix_hybrid' if args.hybrid else 'pourbaix_surface'
    suffix = ''
    if args.gc:
        suffix += '_gc'
    if args.legend_in:
        suffix += '_legend_in'
    elif args.legend_out:
        suffix += '_legend_out'
    elif args.legend_up:
        suffix += '_legend_up'
    if args.suffix:
        suffix += '_' + args.suffix
    return png_name, suffix
    
def load_labels(args, json_files):
    """Load labels and optional corrections from label.csv."""
    label_csv_path = args.label_csv or os.path.join(args.csv_dir, 'label.csv')
    file_labels = {}
    file_gibbs_corrections = {}
    file_gc_params = {}

    if os.path.exists(label_csv_path):
        label_df = pd.read_csv(
            label_csv_path, header=None,
            names=['json_name', 'label', 'G_corr', 'A', 'B', 'C'],
        )
        for _, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']
            if args.gibbs and pd.notna(row.get('G_corr')):
                file_gibbs_corrections[json_name] = float(row['G_corr'])
            if args.gc and all(pd.notna(row.get(col)) for col in ('A', 'B', 'C')):
                file_gc_params[json_name] = {
                    'A': float(row['A']), 'B': float(row['B']), 'C': float(row['C']),
                }
    else:
        for json_file in json_files:
            file_labels[os.path.basename(json_file)] = read(json_file).get_chemical_formula()

    return file_labels, file_gibbs_corrections, file_gc_params

if __name__ == "__main__":
    init_thermo_constants()
    cli_args = parse_args()
    output_png_name, output_suffix = build_output_names(cli_args)
    main(cli_args, output_png_name, output_suffix)