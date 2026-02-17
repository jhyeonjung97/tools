#!/usr/bin/env python3
"""
Example 3: Parse E-fermi from OUTCAR and add it to an existing JSON.

For a directory with both OUTCAR and an ASE JSON file after a VASP run,
read E-fermi from OUTCAR and store it in the JSON row.
"""
import ase.db
import argparse


def get_e_fermi_from_occar(outcar_path='OUTCAR'):
    """Parse the E-fermi line from OUTCAR and return it as a float."""
    with open(outcar_path) as f:
        text = f.read()
    # Use the last occurrence (typically the converged final value)
    return float(text.split('E-fermi :')[-1].split()[0])


def main():
    parser = argparse.ArgumentParser(description='Store E-fermi from OUTCAR into a JSON file')
    parser.add_argument('--json', '-j', default='final_with_calculator.json',
                        help='Path to the target JSON file')
    parser.add_argument('--occar', '-o', default='OUTCAR',
                        help='Path to OUTCAR file')
    parser.add_argument('--id', type=int, default=1,
                        help='Row id in the JSON DB to update')
    args = parser.parse_args()

    e_fermi = get_e_fermi_from_occar(args.occar)
    print(f'OUTCAR E-fermi: {e_fermi} eV')

    with ase.db.connect(args.json, serial=True) as db:
        db.update(args.id, e_fermi=e_fermi)
    print(f'Saved: {args.json} (id={args.id})  e_fermi={e_fermi}')


if __name__ == '__main__':
    main()
