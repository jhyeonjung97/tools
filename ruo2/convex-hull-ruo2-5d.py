from ase.io import read
import os
import pandas as pd
import matplotlib.pyplot as plt
from math import gcd
import mendeleev
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--base-path', type=str, default='~/Desktop/3_RuO2/4_high_valence', help='Base path')
    parser.add_argument('-x', '--fig-width', type=float, default=5.2, help='Figure width (inches)')
    parser.add_argument('-y', '--fig-height', type=float, default=3.6, help='Figure height (inches)')
    parser.add_argument('-o', '--oxygen', type=float, default=-5.29, help='Oxygen chemical potential (eV/atom)')
    parser.add_argument('-t', '--temperature', type=float, default=None, help='Temperature (K)')
    parser.add_argument('--show', action='store_true', help='Show figure')
    args = parser.parse_args()
    base_path = args.base_path
    base_path = os.path.expanduser(base_path)
    fig_width = args.fig_width
    fig_height = args.fig_height
    oxygen = args.oxygen
    temperature = args.temperature
    show = args.show
    if temperature is not None:
        if temperature >= 1000:
            oxygen -= 1.10
        elif temperature >= 900:
            oxygen -= 0.98
        elif temperature >= 800:
            oxygen -= 0.85
        elif temperature >= 700:
            oxygen -= 0.73
        elif temperature >= 600:
            oxygen -= 0.61
        elif temperature >= 500:
            oxygen -= 0.50
        elif temperature >= 400:
            oxygen -= 0.38
        elif temperature >= 300:
            oxygen -= 0.27
        elif temperature >= 200:
            oxygen -= 0.17
        elif temperature >= 100:
            oxygen -= 0.08

    base_path = "~/Desktop/3_RuO2"
    base_path = os.path.expanduser(base_path)
    rows = ['3d', '4d', '5d']
    elements = {
        '3d': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
        '4d': ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'],
        '5d': ['La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    }

    atoms_ruo2 = os.path.join(base_path, '4_high_valence/1_M-RuO2/0_Ru/final_with_calculator.json')
    energy_ruo2 = read(atoms_ruo2).get_potential_energy()
    atoms_reo2 = os.path.join(base_path, '3_ReRuO2_OER/0_bulk/convex_hull/3_ReO2/final_with_calculator.json')
    energy_reo2 = read(atoms_reo2).get_potential_energy()

    for row in rows:
        cmap = plt.colormaps['YlOrRd']
        for i, element in enumerate(elements[row]):
            atomic_number = mendeleev.element(element).atomic_number
            colors = [cmap(i / (len(elements[row]) - 1)) for i in range(len(elements[row]))]
            atoms_mruo2 = os.path.join(base_path, f'4_high_valence/4_M-RuO2_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            atoms_mo2 = os.path.join(base_path, f'4_high_valence/5_MO2_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            atoms_mxoy = os.path.join(base_path, f'4_high_valence/6_MxOy_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            energy_mruo2 = read(atoms_mruo2).get_potential_energy()
            energy_mo2 = read(atoms_mo2).get_potential_energy()
            energy_mxoy = read(atoms_mxoy).get_potential_energy()
            element_count = read(atoms_mxoy).get_chemical_symbols().count(element)
            oxygen_count = read(atoms_mxoy).get_chemical_symbols().count('O')
            element_ratio, oxygen_ratio = get_simplest_formula(element_count, oxygen_count)
            if element_ratio == 1 and oxygen_ratio == 1:
                formula = f'{element}O'
            elif element_ratio == 1:
                formula = f'{element}O$_{oxygen_ratio}$'
            else:
                formula = f'{element}$_{element_ratio}$O$_{oxygen_ratio}$'
            formation_energy_mruo2 = (8*energy_mruo2 - 7*energy_ruo2 - 1*energy_mo2)/8/8
            formation_energy_mxoy = (energy_mxoy/element_count - energy_mo2/8 - oxygen*(oxygen_count/element_count-2))/8
            # if element == 'Co':
            #     print(f"formation_energy_mxoy: {formation_energy_mxoy}")
            #     print(f"energy_mxoy: {energy_mxoy}")
            #     print(f"element_count: {element_count}")
            #     print(f"oxygen_count: {oxygen_count}")
            #     print(f"element_ratio: {element_ratio}")
            #     print(f"oxygen_ratio: {oxygen_ratio}")
            #     print(f"energy_mo2: {energy_mo2}")
            #     print(f"energy_mruo2: {energy_mruo2}")
            #     print(f"energy_ruo2: {energy_ruo2}")
            #     print(f"formation_energy_mxoy: {formation_energy_mxoy}")
            #     print(f"formula: {formula}")
            plt.figure(figsize=(fig_width, fig_height))
            plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green')
            plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green')
            plt.text(0.0, 0.02, 'RuO$_2$\n(rutile)', ha='center', va='bottom', linespacing=0.8)
            plt.text(1.0, 0.02, f'{element}O$_2$\n(rutile)', ha='center', va='bottom', linespacing=0.8)
            plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=0)
            plt.scatter(7/8, formation_energy_mruo2, marker='D', edgecolor='black', facecolor='orange')
            if element == 'Os':
                plt.text(7/8-0.12, formation_energy_mruo2-0.04, f'{element}-RuO$_2$\n(rutile)', ha='center', va='center', linespacing=0.8)
            else:
                plt.text(7/8-0.12, formation_energy_mruo2, f'{element}-RuO$_2$\n(rutile)', ha='center', va='center', linespacing=0.8)
            plt.scatter(1.0, formation_energy_mxoy, marker='s', edgecolor='black', facecolor='blue')
            if element == 'Re':
                formation_energy_reo2 = (energy_reo2/4 - energy_mo2/8)/8
                plt.scatter(1.0, formation_energy_reo2, marker='s', edgecolor='black', facecolor='blue')
                plt.text(1.03, formation_energy_reo2, f'ReO$_2$', ha='left', va='center')
            plt.text(1.03, formation_energy_mxoy, f'{formula}', ha='left', va='center')
            plt.xlabel(r'x, Ru$_x$' + f'{element}' + r'$_{1-x}$O$_y$')
            plt.ylabel('Formation energy (eV/unit)')
            plt.xlim(-0.1, 1.2)
            plt.ylim(-0.4, 0.4)
            plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            plt.savefig(f'RuO2_5d_{atomic_number}_{element}xOy_convex_hull.png', dpi=300, bbox_inches='tight')
            plt.tight_layout()
            if show:
                plt.show()
            plt.close()

def get_simplest_formula(element_count, oxygen_count):
    """원자 개수로부터 가장 간단한 화학식을 구하는 함수
    예: element_count=8, oxygen_count=20 -> (2, 5) -> 'M2O5'
    """
    if element_count == 0 or oxygen_count == 0:
        return None, None
    
    # 최대공약수 구하기
    common_divisor = gcd(element_count, oxygen_count)
    
    # 최대공약수로 나누어 가장 간단한 정수 비 구하기
    element_ratio = element_count // common_divisor
    oxygen_ratio = oxygen_count // common_divisor
    
    return element_ratio, oxygen_ratio

if __name__ == '__main__':
    main()