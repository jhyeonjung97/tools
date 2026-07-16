#!/usr/bin/env python3

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
from ase.io import read

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.sans-serif"] = ["Helvetica"]

ZN_REFERENCE_ENERGY = -1.2363842773  # eV per Zn atom

# MnO2 reference square markers: (phase key after normalize_phase_name) -> offset & alignment
MNO2_ANNOTATIONS = [
    {"phase": "A", "xytext": (-6, 0), "va": "top", "ha": "right"},
    {"phase": "B", "xytext": (-6, 0), "va": "center", "ha": "right"},
    {"phase": "D", "xytext": (-6, 0), "va": "bottom", "ha": "right"},
    {"phase": "E", "xytext": (-6, 0), "va": "center", "ha": "right"},
    {"phase": "G", "xytext": (-6, 0), "va": "bottom", "ha": "right"},
    {"phase": "Spinel", "xytext": (-6, 0), "va": "center", "ha": "right"},
    {"phase": "R", "xytext": (-6, 0), "va": "center", "ha": "right"},
    {"phase": "H", "xytext": (-6, 0), "va": "center", "ha": "right"},
    {"phase": "NO_PHASE", "xytext": (-6, 0), "va": "center", "ha": "right"},
]
MNO2_ANNOTATION_DEFAULT = {"xytext": (-6, 0), "va": "center", "ha": "right"}
MNO2_ANNOTATION_BY_PHASE = {row["phase"]: row for row in MNO2_ANNOTATIONS}

# Selected hull scatter points: JSON stem (lowercase) + text offset / alignment
ZN_STRUCTURE_ANNOTATIONS = [
    {"stem": "e-znmn4o8", "xytext": (-6, 0), "va": "top", "ha": "right"},
    {"stem": "h-znmn2o4", "xytext": (-6, 0), "va": "top", "ha": "right"},
    {"stem": "sp-znmn2o4", "xytext": (-6, 0), "va": "top", "ha": "right"},
    {"stem": "d-znmno2", "xytext": (6, 0), "va": "top", "ha": "left"},
]


def get_counts(atoms):
    return Counter(atoms.get_chemical_symbols())


def normalize_phase_name(phase: str) -> str:
    """Normalize aliases used in filenames."""
    if phase in {"Sp", "Spinel", "H"}:
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
        return "Sp-ZnMn2O4 (Hetaerolite)"

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
        default=str(Path.cwd()),
        help="Directory containing sherry JSON files (default: current working directory)",
    )
    parser.add_argument(
        "--include-no-dash",
        action="store_true",
        help="Include files whose JSON filename has no '-' (default: excluded).",
    )
    parser.add_argument(
        "--palette",
        type=str,
        default="colorblind",
        help="Color palette name (seaborn or matplotlib, default: colorblind)",
    )
    parser.add_argument(
        "--half",
        action="store_true",
        help="Narrow figure width (half) and show x-axis only up to 0.5.",
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

    fig_w = 5 if args.half else 8
    x_max = 0.5 if args.half else 1.1
    plt.figure(figsize=(fig_w, 6))
    plt.axvspan(
        0.33,
        x_max,
        facecolor="whitesmoke",
        zorder=-2,
        linewidth=0,
        edgecolor="none",
    )
    phase_order = ["A", "B", "D", "E", "G", "NO_PHASE", "R", "H", "Sp"]
    phases = sorted({p["phase"] for p in points}, key=lambda x: phase_order.index(x) if x in phase_order else len(phase_order))
    palette = sns.color_palette(args.palette, len(phases))
    phase_colors = {phase: palette[i] for i, phase in enumerate(phases)}
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

        # Draw phase-specific lower convex hull instead of simple line connections.
        phase_hull_candidates = list(phase_points)
        if phase in mno2_by_phase:
            phase_hull_candidates.append(
                {
                    "name": f"{phase}-MnO2_ref",
                    "x": 0.0,
                    "y": mno2_by_phase[phase] - mn02_ref_energy,
                }
            )
        phase_hull_candidates.append({"name": "Zn_ref", "x": 1.0, "y": 0.0})
        phase_hull = lower_hull(phase_hull_candidates)
        if len(phase_hull) >= 2:
            plt.plot(
                [p["x"] for p in phase_hull],
                [p["y"] for p in phase_hull],
                color=phase_colors[phase],
                lw=1.2,
                alpha=0.9,
                zorder=2,
            )

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
            mn02_label = f"{display_phase_name(phase)}-MnO2"
            style = MNO2_ANNOTATION_BY_PHASE.get(phase, MNO2_ANNOTATION_DEFAULT)
            plt.annotate(
                subscript_digits(mn02_label),
                (0.0, y_mno2),
                xytext=style["xytext"],
                textcoords="offset points",
                fontsize=10,
                alpha=0.95,
                va=style["va"],
                ha=style["ha"],
            )

    zn_structure_annotations = (
        [
            {"stem": "e-znmn4o8", "xytext": (-6, 0), "va": "top", "ha": "right"},
            {"stem": "h-znmn2o4", "xytext": (-6, 0), "va": "top", "ha": "right"},
            {"stem": "sp-znmn2o4", "xytext": (-6, 0), "va": "top", "ha": "right"},
        ]
        if args.half
        else ZN_STRUCTURE_ANNOTATIONS
    )
    zn_ann_by_stem = {row["stem"].lower(): row for row in zn_structure_annotations}
    for p in points:
        row = zn_ann_by_stem.get(p["name"].lower())
        if row is None:
            continue
        plt.annotate(
            subscript_digits(display_structure_name(p["name"])),
            (p["x"], p["y"]),
            xytext=row["xytext"],
            textcoords="offset points",
            fontsize=10,
            alpha=0.9,
            va=row["va"],
            ha=row["ha"],
        )

    hx = [p["x"] for p in hull_points]
    hy = [p["y"] for p in hull_points]
    plt.plot(hx, hy, "--", color="black", lw=1.0, label="overall", zorder=0)
    plt.scatter([1.0], [0.0], s=55, marker="s", facecolors="silver", edgecolors="black", linewidths=0.3, zorder=4)
    plt.annotate(
        subscript_digits("Zn"),
        (1.0, 0.0),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=10,
        alpha=0.95,
        va="center",
    )

    plt.plot([0.0, 1.0], [0.0, 0.0], color="silver", lw=0.5, linestyle="-", zorder=-1)
    plt.xlim(-0.1, x_max)
    if args.half:
        plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    plt.ylim(-0.9, 0.3)
    plt.tick_params(labelsize=10)
    plt.xlabel(r"Zn composition, Zn$_x$(MnO$_2$)$_{1-x}$", fontsize=12)
    plt.ylabel(r"Formation energy ($\Delta E$, eV per $n_\mathrm{Mn}+n_{\mathrm{Zn}}$)", fontsize=12)
    plt.legend(loc="upper right" if args.half else "lower right", fontsize=10)
    plt.tight_layout()
    output_path = sherry_dir / f"convex-hull-mno2-zn-{args.palette}{'-half' if args.half else ''}.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
