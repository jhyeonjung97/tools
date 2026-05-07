#!/usr/bin/env python3

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from ase.io import read


ZN_REFERENCE_ENERGY = -1.2363842773  # eV per Zn atom
ANNOTATION_INCLUDE_ZN_STEMS = {
    "sp-znmn4o8",
    "e-znmn4o8",
    "e-znmn2o4",
    "r-znmn2o4",
    "h-znmn2o4",
    "sp-znmn2o4",
    "sp-znmno2",
    "d-znmno2",
}


def get_counts(atoms):
    return Counter(atoms.get_chemical_symbols())


def normalize_phase_name(phase: str) -> str:
    """Normalize aliases used in filenames."""
    if phase in {"Sp", "Spinel"}:
        return "Spinel"
    return phase


def display_phase_name(phase: str) -> str:
    """Display label for phase names in legend/annotations."""
    phase_map = {
        "A": "α",
        "B": "β",
        "D": "δ",
        "E": "ε",
        "Spinel": "λ",
        "G": "ɣ",
    }
    return phase_map.get(phase, phase)


def display_structure_name(name: str) -> str:
    """Display label for structure names in annotations."""
    if name.lower() == "h-znmn2o4":
        return "Hetaerolite (ZnMn2O4)"

    prefix_map = {
        "A-": "α-",
        "B-": "β-",
        "D-": "δ-",
        "E-": "ε-",
        "Sp-": "λ-",
        "Spinel-": "λ-",
        "G-": "ɣ-",
    }
    for src, dst in prefix_map.items():
        if name.startswith(src):
            return dst + name[len(src):]
    return name


def subscript_digits(text: str) -> str:
    """Convert digits to Unicode subscripts for annotation labels."""
    table = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return text.translate(table)


def reduced_formula_key(counts):
    """Return reduced Zn-Mn-O composition key, e.g. Zn1Mn2O4."""
    n_zn = counts.get("Zn", 0)
    n_mn = counts.get("Mn", 0)
    n_o = counts.get("O", 0)
    if n_zn == 0 and n_mn == 0:
        return "INVALID"

    from math import gcd

    g = gcd(gcd(max(n_zn, 0), max(n_mn, 0)), max(n_o, 0))
    if g == 0:
        g = 1
    r_zn, r_mn, r_o = n_zn // g, n_mn // g, n_o // g

    return f"Zn{r_zn}Mn{r_mn}O{r_o}"


def find_mno2_reference_energy(sherry_dir: Path):
    """Find the most stable MnO2 reference among JSON files in sherry_dir."""
    mno2_candidates = []
    for path in sorted(sherry_dir.glob("*.json")):
        atoms = read(path)
        counts = get_counts(atoms)

        # Exact MnO2 stoichiometry only (no Zn/H/etc.)
        if set(counts.keys()) != {"Mn", "O"}:
            continue
        n_mn = counts["Mn"]
        n_o = counts["O"]
        if n_mn <= 0 or n_o != 2 * n_mn:
            continue

        energy_per_mno2 = atoms.get_potential_energy() / n_mn
        mno2_candidates.append((energy_per_mno2, path.name))

    if not mno2_candidates:
        raise RuntimeError("No exact MnO2 structures found in sherry directory.")

    mno2_candidates.sort(key=lambda x: x[0])
    return mno2_candidates[0], mno2_candidates


def build_points(sherry_dir: Path, mn02_ref_energy: float, include_no_dash: bool = False):
    """Build (x_Zn, formation_energy) points from *Zn*.json files."""
    excluded_files = {"C-ZnMn3O7.json"}
    points = []
    for path in sorted(sherry_dir.glob("*Zn*.json")):
        if path.name in excluded_files:
            continue
        if (not include_no_dash) and ("-" not in path.stem):
            continue

        atoms = read(path)
        counts = get_counts(atoms)

        n_zn = counts.get("Zn", 0)
        n_mn = counts.get("Mn", 0)
        if n_zn <= 0 or n_mn <= 0:
            continue

        total_metal = n_zn + n_mn
        x_zn = n_zn / total_metal
        total_energy = atoms.get_potential_energy()

        # Reference: n_Mn * E(MnO2_min) + n_Zn * E(Zn_metal_ref)
        formation_energy = (
            total_energy - n_mn * mn02_ref_energy - n_zn * ZN_REFERENCE_ENERGY
        ) / total_metal

        points.append(
            {
                "name": path.stem,
                "path": str(path),
                "phase": normalize_phase_name(path.stem.split("-", 1)[0]) if "-" in path.stem else "NO_PHASE",
                "composition_key": reduced_formula_key(counts),
                "x": x_zn,
                "y": formation_energy,
            }
        )

    if not points:
        raise RuntimeError("No valid Zn-containing Mn structures found in sherry directory.")
    return points


def lower_hull(points):
    """Return lower convex hull points from [{'x','y',...}, ...]."""
    by_x = {}
    for p in points:
        x = p["x"]
        if x not in by_x or p["y"] < by_x[x]["y"]:
            by_x[x] = p

    sorted_points = sorted(by_x.values(), key=lambda p: p["x"])

    def cross(a, b, c):
        return (b["x"] - a["x"]) * (c["y"] - a["y"]) - (b["y"] - a["y"]) * (c["x"] - a["x"])

    hull = []
    for p in sorted_points:
        while len(hull) >= 2 and cross(hull[-2], hull[-1], p) <= 0:
            hull.pop()
        hull.append(p)
    return hull


def main():
    parser = argparse.ArgumentParser(
        description="Build Zn-composition convex hull from mno2 sherry JSON files."
    )
    parser.add_argument(
        "--sherry-dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "figures" / "sherry"),
        help="Directory containing sherry JSON files (default: ./figures/sherry)",
    )
    parser.add_argument(
        "--include-no-dash",
        action="store_true",
        help="Include files whose JSON filename has no '-' (default: excluded).",
    )
    args = parser.parse_args()

    sherry_dir = Path(args.sherry_dir).expanduser().resolve()
    if not sherry_dir.exists():
        raise FileNotFoundError(f"Directory not found: {sherry_dir}")

    (mn02_ref_energy, mn02_ref_name), mn02_candidates = find_mno2_reference_energy(sherry_dir)
    points = build_points(sherry_dir, mn02_ref_energy, include_no_dash=args.include_no_dash)

    # Add synthetic reference endpoints for visual hull context.
    endpoint_mno2 = {"name": f"MnO2_ref ({mn02_ref_name})", "x": 0.0, "y": 0.0}
    endpoint_zn = {"name": "Zn_ref", "x": 1.0, "y": 0.0}
    hull_points = lower_hull(points + [endpoint_mno2, endpoint_zn])

    print(f"Using sherry directory: {sherry_dir}")
    print(f"Zn reference energy: {ZN_REFERENCE_ENERGY:.10f} eV/Zn")
    print(f"MnO2 reference (most stable): {mn02_ref_name}")
    print(f"MnO2 reference energy: {mn02_ref_energy:.10f} eV per MnO2")
    print(f"Loaded Zn-containing structures: {len(points)}")
    print("MnO2 candidates (eV per MnO2):")
    for energy, name in mn02_candidates:
        print(f"  {name:30s} {energy: .8f}")

    plt.figure(figsize=(10, 6))
    phases = sorted({p["phase"] for p in points})
    cmap = plt.get_cmap("tab20")
    phase_colors = {phase: cmap(i % 20) for i, phase in enumerate(phases)}
    mno2_by_phase = {}
    for energy, name in mn02_candidates:
        stem = name.replace(".json", "")
        if "-" in stem:
            phase = normalize_phase_name(stem.split("-", 1)[0])
            mno2_by_phase[phase] = energy

    for phase in phases:
        phase_points = [p for p in points if p["phase"] == phase]
        phase_points = sorted(phase_points, key=lambda p: p["x"])
        xs = [p["x"] for p in phase_points]
        ys = [p["y"] for p in phase_points]

        # Connect phase-specific MnO2 -> phase points -> shared Zn endpoint.
        # Skip line connection for H phase.
        if phase != "H":
            line_xs = []
            line_ys = []
            if phase in mno2_by_phase:
                line_xs.append(0.0)
                line_ys.append(mno2_by_phase[phase] - mn02_ref_energy)
            line_xs.extend(xs)
            line_ys.extend(ys)
            line_xs.append(1.0)
            line_ys.append(0.0)
            plt.plot(line_xs, line_ys, color=phase_colors[phase], lw=1.0, alpha=0.85, zorder=2)

        plt.scatter(
            xs,
            ys,
            s=45,
            color=phase_colors[phase],
            edgecolors="black",
            linewidths=0.3,
            label=f"{display_phase_name(phase)} phase" if phase != "H" else "_nolegend_",
            zorder=3,
        )

        # Mark each phase's MnO2 point explicitly.
        if phase in mno2_by_phase:
            y_mno2 = mno2_by_phase[phase] - mn02_ref_energy
            plt.scatter(
                [0.0],
                [y_mno2],
                s=50,
                marker="s",
                color=phase_colors[phase],
                edgecolors="black",
                linewidths=0.3,
                zorder=4,
            )
            plt.annotate(
                subscript_digits(f"{display_phase_name(phase)}-MnO2"),
                (0.0, y_mno2),
                xytext=(-6, 0) if phase in {"G", "R"} else (6, 0),
                textcoords="offset points",
                fontsize=9,
                alpha=0.95,
                va="center",
                ha="right" if phase in {"G", "R"} else "left",
            )

    for p in points:
        if p["name"].lower() not in ANNOTATION_INCLUDE_ZN_STEMS:
            continue
        is_left_label = p["name"].lower() == "sp-znmno2"
        label_offset = (-6, 0) if is_left_label else (6, 0)
        plt.annotate(
            subscript_digits(display_structure_name(p["name"])),
            (p["x"], p["y"]),
            xytext=label_offset,
            textcoords="offset points",
            fontsize=9,
            alpha=0.9,
            va="center",
            ha="right" if is_left_label else "left",
        )

    hx = [p["x"] for p in hull_points]
    hy = [p["y"] for p in hull_points]
    plt.plot(hx, hy, "--", color="black", lw=1.0, label="lower convex hull", zorder=0)
    plt.scatter([1.0], [0.0], s=55, marker="s", facecolors="white", edgecolors="black", linewidths=0.3, zorder=4)
    plt.annotate(
        subscript_digits("Zn"),
        (1.0, 0.0),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=9,
        alpha=0.95,
        va="center",
    )

    plt.axhline(0.0, color="silver", lw=0.8, ls="--", zorder=-1)
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.9, 0.3)
    plt.xlabel("Zn composition, Zn$_x$Mn$_{1-x}$O$_2$")
    plt.ylabel("Formation energy (ΔE, eV/metal)")
    plt.legend()
    plt.tight_layout()
    output_path = sherry_dir / "convex-hull-mno2-zn.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
