#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read

ANION_SYMBOLS = {"O", "H", "F", "Cl", "Br", "I", "N", "S"}


def metal_symbols(atoms) -> list[str]:
    """Return metal element symbols present in the structure."""
    return sorted({s for s in atoms.get_chemical_symbols() if s not in ANION_SYMBOLS})


def mo_distances(atoms, metal: str | None = None, cutoff: float = 2.5) -> list[float]:
    """Return all M-O distances within cutoff (Å), using periodic boundary conditions."""
    symbols = atoms.get_chemical_symbols()
    o_indices = [i for i, s in enumerate(symbols) if s == "O"]
    if not o_indices:
        raise ValueError("No oxygen atoms found in structure.")

    metals = [metal] if metal else metal_symbols(atoms)
    if not metals:
        raise ValueError("No metal atoms found in structure.")

    distances = []
    for m in metals:
        m_indices = [i for i, s in enumerate(symbols) if s == m]
        if not m_indices:
            if metal:
                raise ValueError(f"Metal '{metal}' not found in structure.")
            continue

        for i in m_indices:
            for j in o_indices:
                d = atoms.get_distance(i, j, mic=True)
                if d <= cutoff:
                    distances.append(d)

    if not distances:
        raise RuntimeError(
            f"No M-O pairs found within {cutoff:.2f} Å. Try increasing --cutoff."
        )

    return distances


def average_mo_bond_length(atoms, metal: str | None = None, cutoff: float = 2.5) -> dict:
    """Compute average M-O bond length using pairs within cutoff (Å)."""
    distances = mo_distances(atoms, metal=metal, cutoff=cutoff)
    arr = np.array(distances)
    return {
        "n_bonds": len(arr),
        "mean": float(arr.mean()),
        "std": float(arr.std()),
        "min": float(arr.min()),
        "max": float(arr.max()),
    }


def mean_mo_bond_length(atoms, metal: str | None = None, cutoff: float = 2.5) -> float:
    """Return pooled mean M-O bond length (Å)."""
    return average_mo_bond_length(atoms, metal=metal, cutoff=cutoff)["mean"]


def main():
    parser = argparse.ArgumentParser(
        description="Compute average M-O bond length from ASE JSON structure files."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="JSON file or directory (default: current directory)",
    )
    parser.add_argument(
        "--metal",
        type=str,
        default=None,
        help="Metal symbol (default: auto-detect all non-anion elements)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=2.5,
        help="Maximum M-O distance counted as a bond, in Å (default: 2.5)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed statistics for each structure",
    )
    args = parser.parse_args()

    target = Path(args.path).expanduser().resolve()
    if target.is_file():
        json_paths = [target]
    elif target.is_dir():
        json_paths = sorted(target.glob("*.json"))
        if not json_paths:
            print(f"No JSON files found in: {target}", file=sys.stderr)
            sys.exit(1)
    else:
        raise FileNotFoundError(f"Path not found: {target}")

    for json_path in json_paths:
        try:
            atoms = read(json_path)
            stats = average_mo_bond_length(atoms, metal=args.metal, cutoff=args.cutoff)
        except Exception as exc:
            print(f"{json_path.stem} ERROR: {exc}", file=sys.stderr)
            continue

        if args.verbose:
            print(f"Structure: {json_path.stem}")
            print(f"  mean = {stats['mean']:.4f} Å")
            print(f"  std  = {stats['std']:.4f} Å")
            print(f"  min  = {stats['min']:.4f} Å")
            print(f"  max  = {stats['max']:.4f} Å")
            print(f"  n    = {stats['n_bonds']}")
        else:
            print(f"{json_path.stem} {stats['mean']:.4f}")


if __name__ == "__main__":
    main()
