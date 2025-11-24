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

    base_path = "~/Desktop/3_RuO2/4_high_valence"
    base_path = os.path.expanduser(base_path)
    rows = ['3d', '4d', '5d']
    elements = {
        '3d': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
        '4d': ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'],
        '5d': ['La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    }

    atoms_ruo2 = os.path.join(base_path, '1_M-RuO2/0_Ru/final_with_calculator.json')
    energy_ruo2 = read(atoms_ruo2).get_potential_energy()

    for row in rows:
        plt.figure(figsize=(fig_width, fig_height))
        plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green')
        plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green')
        plt.text(0.0, 0.01, 'RuO$_2$', ha='center', va='bottom')
        plt.text(1.0, 0.01, 'MO$_2$', ha='center', va='bottom')
        plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=0)
        cmap = plt.colormaps['YlOrRd']
        for i, element in enumerate(elements[row]):
            atomic_number = mendeleev.element(element).atomic_number
            colors = [cmap(i / (len(elements[row]) - 1)) for i in range(len(elements[row]))]
            atoms_mo2 = os.path.join(base_path, f'5_MO2_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            atoms_mruo2 = os.path.join(base_path, f'4_M-RuO2_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            energy_mo2 = read(atoms_mo2).get_potential_energy()
            energy_mruo2 = read(atoms_mruo2).get_potential_energy()
            formation_energy = (8*energy_mruo2 - 7*energy_ruo2 - 1*energy_mo2)/8/8
            if element != 'La':
                plt.scatter(7/8, formation_energy, marker='D', edgecolor='black', facecolor=colors[i])
            if element in ['Cu', 'Y', 'Mo', 'Re', 'Pt']:
                plt.text(7/8+0.03, formation_energy, f'{element}-RuO$_2$', ha='left', va='center')
            elif element in ['Tc', 'Ir']:
                plt.text(7/8-0.03, formation_energy+0.01, f'{element}-RuO$_2$', ha='right', va='center')
            elif element in ['Ru', 'Nb']:
                plt.text(7/8-0.03, formation_energy-0.01, f'{element}-RuO$_2$', ha='right', va='center')
            elif element != 'La':
                plt.text(7/8-0.03, formation_energy, f'{element}-RuO$_2$', ha='right', va='center')
        plt.xlabel(r'x, Ru$_x$M$_{1-x}$O$_2$')
        plt.ylabel('Formation energy (eV/unit)')
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.3, 0.3)
        plt.savefig(f'RuO2_{row}_MO2_convex_hull.png', dpi=300, bbox_inches='tight')
        plt.tight_layout()
        if show:
            plt.show()
        plt.close()

    for row in rows:
        for i, element in enumerate(elements[row]):
            atomic_number = mendeleev.element(element).atomic_number
            atoms_mruo2 = os.path.join(base_path, f'4_M-RuO2_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            atoms_mxoy = os.path.join(base_path, f'6_MxOy_all/precise2/{atomic_number}_{element}/final_with_calculator.json')
            if os.path.exists(atoms_mruo2) and os.path.exists(atoms_mxoy):
                energy_mruo2 = read(atoms_mruo2).get_potential_energy()
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
                formation_energy = (energy_mruo2 - energy_ruo2*7/8 - energy_mxoy/element_count - oxygen*(2-oxygen_count/element_count))/8
                plt.figure(figsize=(fig_width, fig_height))
                plt.scatter(0.0, 0.0, marker='s', edgecolor='black', facecolor='green')
                plt.scatter(1.0, 0.0, marker='s', edgecolor='black', facecolor='green')
                plt.text(0.0, 0.01, 'RuO$_2$', ha='center', va='bottom')
                plt.text(1.0, 0.01, formula, ha='center', va='bottom')
                plt.plot([0.0, 1.0], [0.0, 0.0], color='black', linestyle='-', zorder=0)
                plt.scatter(7/8, formation_energy, marker='D', edgecolor='black', facecolor='orange')
                plt.text(7/8-0.03, formation_energy, f'{element}-RuO$_2$', ha='right', va='center')
                plt.xlabel(r'x, Ru$_x$M$_{1-x}$O$_y$')
                plt.ylabel('Formation energy (eV/unit)')
                plt.xlim(-0.1, 1.1)
                plt.ylim(-0.3, 0.3)
                plt.savefig(f'RuO2_{atomic_number}_{element}xOy_convex_hull.png', dpi=300, bbox_inches='tight')
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