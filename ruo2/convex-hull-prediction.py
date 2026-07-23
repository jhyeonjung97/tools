from ase.io import read
import os
import matplotlib.pyplot as plt
from math import gcd
import mendeleev
import argparse

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

# dH of molecular references at 300K:
# 1/2 O2 derived from H2O/H2:  -4.658724749999999
# 1/2 O2 derived from PBE:  -4.851213114999999

# Host definitions: name -> (M-host subdir, reference element for pure host oxide)
HOSTS = {
    'RuO2': {'dir': '3_M-RuO2', 'ref': '44_Ru', 'label': 'RuO$_2$', 'dope': 'RuO$_2$'},
    'IrO2': {'dir': '4_M-IrO2', 'ref': '77_Ir', 'label': 'IrO$_2$', 'dope': 'IrO$_2$'},
}

elements = {
    '3d': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    '4d': ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'],
    '5d': ['La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--base-path', type=str,
                        default='~/Desktop/7_prediction/1_bulk_opt',
                        help='Base path of bulk optimization data')
    parser.add_argument('-d', '--output-dir', type=str,
                        default='~/Desktop/7_prediction/figures/convex_hulls',
                        help='Output directory for figures')
    parser.add_argument('--host', type=str, default='both',
                        choices=['RuO2', 'IrO2', 'both'],
                        help='Which host oxide convex hull(s) to plot')
    parser.add_argument('-x', '--fig-width', type=float, default=3.2, help='Figure width (inches)')
    parser.add_argument('-y', '--fig-height', type=float, default=3.2, help='Figure height (inches)')
    parser.add_argument('-o', '--oxygen', type=float,
                        default=-4.658724749999999 + 0.27 - 0.73,
                        help='Oxygen chemical potential (eV/atom)')
    parser.add_argument('-t', '--temperature', type=float, default=None, help='Temperature (K)')
    parser.add_argument('--show', action='store_true', help='Show figure')
    args = parser.parse_args()

    base_path = os.path.expanduser(args.base_path)
    output_dir = os.path.expanduser(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    oxygen = args.oxygen
    temperature = args.temperature
    if temperature is not None:
        # temperature-dependent oxygen entropy correction (delta mu_O)
        for thr, corr in [(1000, 1.10), (900, 0.98), (800, 0.85), (700, 0.73),
                          (600, 0.61), (500, 0.50), (400, 0.38), (300, 0.27),
                          (200, 0.17), (100, 0.08)]:
            if temperature >= thr:
                oxygen -= corr
                break
    print('oxygen: ', oxygen)

    hosts = ['RuO2', 'IrO2'] if args.host == 'both' else [args.host]
    for host in hosts:
        plot_host(host, base_path, output_dir, oxygen,
                  args.fig_width, args.fig_height, args.show)


def plot_host(host, base_path, output_dir, oxygen, fig_width, fig_height, show):
    info = HOSTS[host]
    host_dir = info['dir']
    print(f'\n=== Host: {host} ===')

    # pure host oxide reference (e.g. RuO2 = Ru in M-RuO2, 8 Ru + 16 O)
    atoms_host = os.path.join(base_path, host_dir, info['ref'], 'final_with_calculator.json')
    energy_host = read(atoms_host).get_potential_energy()

    for row in ['3d', '4d', '5d']:
        for element in elements[row]:
            atomic_number = mendeleev.element(element).atomic_number
            key = f'{atomic_number}_{element}'
            f_mdoped = os.path.join(base_path, host_dir, key, 'final_with_calculator.json')
            f_mo2 = os.path.join(base_path, '2_MO2', key, 'final_with_calculator.json')
            f_mxoy = os.path.join(base_path, '1_MxOy', key, 'final_with_calculator.json')
            if not (os.path.exists(f_mdoped) and os.path.exists(f_mo2) and os.path.exists(f_mxoy)):
                print(f'  skip {key}: missing data')
                continue

            energy_mdoped = read(f_mdoped).get_potential_energy()
            energy_mo2 = read(f_mo2).get_potential_energy()
            atoms_mxoy = read(f_mxoy)
            energy_mxoy = atoms_mxoy.get_potential_energy()
            symbols = atoms_mxoy.get_chemical_symbols()
            element_count = symbols.count(element)
            oxygen_count = symbols.count('O')
            element_ratio, oxygen_ratio = get_simplest_formula(element_count, oxygen_count)
            if element_ratio is None:
                print(f'  skip {key}: zero count')
                continue
            ratio = element_count / element_ratio

            if element_ratio == 1 and oxygen_ratio == 1:
                formula = f'{element}O'
            elif element_ratio == 1:
                formula = f'{element}O$_{{{oxygen_ratio}}}$'
            else:
                formula = f'{element}$_{{{element_ratio}}}$O$_{{{oxygen_ratio}}}$'

            # formation energy of M-doped host (1/8 substitution), per atom
            formation_energy_mdoped = (8 * energy_mdoped - 7 * energy_host - 1 * energy_mo2) / 8 / 24
            # formation energy of MxOy relative to MO2 + O reservoir, per atom
            formation_energy_mxoy = (energy_mxoy / ratio - energy_mo2 / 8 * element_ratio
                                     - oxygen * (oxygen_ratio - 2 * element_ratio)) / element_ratio / 3

            print(f'  {element}: dE(MxOy - M-doped) = {formation_energy_mxoy - formation_energy_mdoped:.3f}')

            plt.figure(figsize=(fig_width, fig_height))
            # endpoints: host oxide (x=0) and MO2 (x=1), both rutile, on the 0 line
            plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green', zorder=2)
            plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green', zorder=2)
            plt.text(0.0, 0.02, f"{info['label']}\n(rutile)", ha='center', va='bottom', linespacing=0.8)
            plt.text(1.0, 0.02, f'{element}O$_2$\n(rutile)', ha='center', va='bottom', linespacing=0.8)
            plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=1)
            # M-doped host at x = 1/8
            plt.scatter(1 / 8, formation_energy_mdoped, marker='D', edgecolor='black', facecolor='orange', zorder=2)
            plt.text(1 / 8, 0.02, f"{element}-{info['dope']}\n(rutile)", ha='left', va='bottom', linespacing=0.8)
            # MxOy at x = 1
            plt.scatter(1.0, formation_energy_mxoy, marker='s', edgecolor='black', facecolor='blue', zorder=2)
            plt.plot([0.0, 1.0], [0.0, formation_energy_mxoy], color='black', linestyle='-', zorder=1)
            plt.plot([1 / 8, 1 / 8], [formation_energy_mdoped, formation_energy_mxoy / 8],
                     color='red', linestyle='-', zorder=0)
            delta = -formation_energy_mxoy / 8 + formation_energy_mdoped
            plt.text(1 / 8 + 0.02, formation_energy_mxoy / 16 + formation_energy_mdoped / 2 - 0.02,
                     f'{delta:+.2f} eV',
                     ha='left', va='center', color='red')
            if formation_energy_mxoy > 0:
                plt.text(1.0, formation_energy_mxoy + 0.03, f'{formula}', ha='center', va='bottom')
            else:
                plt.text(1.0, formation_energy_mxoy - 0.03, f'{formula}', ha='center', va='top')

            host_m = host[:-2]  # 'Ru' or 'Ir'
            plt.xlabel(r'x, ' + f'{element}' + r'$_{x}$' + f'{host_m}' + r'$_{1-x}$O$_y$', fontsize=12)
            plt.ylabel('Formation energy (eV/atom)', fontsize=12)
            plt.xlim(-0.15, 1.15)
            plt.ylim(-1.0, 0.4)
            plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            plt.tight_layout()
            out = os.path.join(output_dir, f'{host}_{atomic_number}_{element}xOy_convex_hull')
            plt.savefig(out + '.png', dpi=300, bbox_inches='tight')
            plt.savefig(out + '.pdf', dpi=300, bbox_inches='tight')
            if show:
                plt.show()
            plt.close()


def get_simplest_formula(element_count, oxygen_count):
    """Return the simplest integer ratio of a MxOy formula.
    e.g. element_count=8, oxygen_count=20 -> (2, 5) -> 'M2O5'
    """
    if element_count == 0 or oxygen_count == 0:
        return None, None
    common_divisor = gcd(element_count, oxygen_count)
    return element_count // common_divisor, oxygen_count // common_divisor


if __name__ == '__main__':
    main()
