#!/usr/bin/env python3
"""LOBSTER ICOHPLIST.lobster 집계.

기본: COHP#(같은 결합 블록)마다
  - 금속 atomMU: s·d 오비탈 행만
  - 산소 atomNU: p 오비탈 행만 (O*_2p_* 등; O*_2s 제외)
  위 조건을 동시에 만족하는 행들의 ICOHP만 합산 (예: COHP#1의 line 4–6, 8–26).

옵션으로 라벨별·site별 예전 집계도 가능.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

ATOM_MU_RE = re.compile(r"^([A-Z][a-z]?)(\d+)(?:_(.+))?$")
METAL_S = re.compile(r"_(\d+)s(?:$|[^0-9])")
METAL_D = re.compile(r"_(\d+)d")
SHELL_P = re.compile(r"_(\d+)p_")


def iter_icohp_rows(path: Path, spin_filter: int | None):
    """Yields: cohp_id, atom_mu, atom_nu, distance, icohp."""
    current_spin = 1
    spin_header = re.compile(r"for spin\s+(\d+)\s*$")

    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            m = spin_header.search(line)
            if m and line.lstrip().startswith("COHP#"):
                current_spin = int(m.group(1))
                continue
            if spin_filter is not None and current_spin != spin_filter:
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            if not parts[0].isdigit():
                continue
            try:
                cohp_id = int(parts[0])
                dist = float(parts[3])
                icohp = float(parts[7])
            except ValueError:
                continue
            yield cohp_id, parts[1], parts[2], dist, icohp


def is_coarse_atom(label: str) -> bool:
    """Ru1, O19 처럼 오비탈 접미사 없는 라벨."""
    m = ATOM_MU_RE.match(label)
    return m is not None and m.group(3) is None


def is_metal_sd_atom_mu(atom_mu: str) -> bool:
    m = ATOM_MU_RE.match(atom_mu)
    if not m or m.group(3) is None:
        return False
    elem = m.group(1)
    if elem == "O":
        return False
    if SHELL_P.search(atom_mu):
        return False
    if re.search(r"_(\d+)f", atom_mu):
        return False
    return bool(METAL_S.search(atom_mu) or METAL_D.search(atom_mu))


def is_oxygen_p_label(atom: str) -> bool:
    m = ATOM_MU_RE.match(atom)
    if not m or m.group(1) != "O" or m.group(3) is None:
        return False
    return bool(SHELL_P.search(atom))


def row_metal_sd_and_oxygen_p(atom_mu: str, atom_nu: str) -> bool:
    """금속 s,d (MU) & 산소 p (NU), 또는 반대 배치."""
    if is_metal_sd_atom_mu(atom_mu) and is_oxygen_p_label(atom_nu):
        return True
    if is_metal_sd_atom_mu(atom_nu) and is_oxygen_p_label(atom_mu):
        return True
    return False


def oxygen_site_key(atom: str) -> str | None:
    m = re.match(r"^(O\d+)_(\d+)p_", atom)
    return m.group(1) if m else None


def oxygen_site_from_row(atom_mu: str, atom_nu: str) -> str | None:
    if is_oxygen_p_label(atom_nu):
        return oxygen_site_key(atom_nu)
    if is_oxygen_p_label(atom_mu):
        return oxygen_site_key(atom_mu)
    return None


def sum_by_cohp_sd_mu_p_nu(
    path: Path, spin_filter: int | None
) -> tuple[dict[int, float], dict[int, tuple[str, str, float]]]:
    """COHP# → (sd 금속 & p 산소) 행만 합산. 메타는 coarse Ru1–O19 행에서."""
    sums: dict[int, float] = defaultdict(float)
    meta: dict[int, tuple[str, str, float]] = {}

    for cohp_id, am, an, dist, icohp in iter_icohp_rows(path, spin_filter):
        if is_coarse_atom(am) and is_coarse_atom(an):
            m_mu = ATOM_MU_RE.match(am)
            m_nu = ATOM_MU_RE.match(an)
            if (
                m_mu
                and m_nu
                and m_mu.group(1) != "O"
                and m_nu.group(1) == "O"
                and cohp_id not in meta
            ):
                meta[cohp_id] = (am, an, dist)
            if (
                m_mu
                and m_nu
                and m_nu.group(1) != "O"
                and m_mu.group(1) == "O"
                and cohp_id not in meta
            ):
                meta[cohp_id] = (am, an, dist)

        if row_metal_sd_and_oxygen_p(am, an):
            sums[cohp_id] += icohp

    return dict(sums), meta


def sum_metal_sd_by_atom_mu(path: Path, spin_filter: int | None) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    for _id, am, _an, _dist, icohp in iter_icohp_rows(path, spin_filter):
        if is_metal_sd_atom_mu(am):
            sums[am] += icohp
    return dict(sums)


def sum_oxygen_p_by_site(path: Path, spin_filter: int | None) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    for _id, am, an, _dist, icohp in iter_icohp_rows(path, spin_filter):
        site = oxygen_site_from_row(am, an)
        if site:
            sums[site] += icohp
    return dict(sums)


def coarse_label_ele_idx(label: str) -> tuple[str, str]:
    """coarse 라벨 Ru1, Hf16 → ('Ru','1'), ('Hf','16'). 실패 시 ('?','?')."""
    m = ATOM_MU_RE.match(label)
    if not m or m.group(3) is not None:
        return ("?", "?")
    return m.group(1), m.group(2)


def sort_key_mu(s: str) -> tuple:
    m = re.match(r"^([A-Za-z]+)(\d+)(?:_|$)", s)
    if m:
        return (m.group(1), int(m.group(2)), s)
    return ("", 0, s)


def sort_key_o(site: str) -> tuple:
    m = re.match(r"^O(\d+)$", site)
    if m:
        return (int(m.group(1)),)
    return (0,)


def main() -> None:
    p = argparse.ArgumentParser(
        description="ICOHPLIST: COHP#별(sd금속×p산소) 합, 기타 옵션"
    )
    p.add_argument(
        "icohplist",
        type=Path,
        nargs="?",
        default="ICOHPLIST.lobster",
        help="ICOHPLIST.lobster 경로",
    )
    p.add_argument(
        "--spin",
        type=int,
        choices=(1, 2),
        default=None,
        help="해당 스핀만 (미지정 시 스핀 1+2)",
    )
    p.add_argument(
        "--by-orbital-label",
        action="store_true",
        help="금속 atomMU(s,d) 라벨별 전역 합 (이전 방식)",
    )
    p.add_argument(
        "--by-oxygen-site",
        action="store_true",
        help="산소 O* site별 p 오비탈 합만 추가 출력",
    )
    p.add_argument(
        "--all-orbitals",
        action="store_true",
        help="오비탈 무시하고 atomMU 문자열만 전부 합산",
    )
    args = p.parse_args()

    if not args.icohplist.is_file():
        raise SystemExit(f"파일 없음: {args.icohplist}")

    spin_note = f" (spin {args.spin}만)" if args.spin else " (spin 1+2)"

    if args.all_orbitals:
        sums: dict[str, float] = defaultdict(float)
        for _id, am, _n, _dist, icohp in iter_icohp_rows(args.icohplist, args.spin):
            sums[am] += icohp
        print(f"# {args.icohplist}{spin_note}  —  atomMU별 sum(ICOH), 전 오비탈")
        print(f"# {'atomMU':<22}  sum_ICOHP")
        for k in sorted(sums, key=sort_key_mu):
            print(f"  {k:<22}  {sums[k]:.5f}")
        return

    sums, meta = sum_by_cohp_sd_mu_p_nu(args.icohplist, args.spin)

    # print(f"# {args.icohplist}{spin_note}")
    # print(
    #     "# COHP#별 합: 금속 s·d (atomMU 또는 NU) × 산소 p (반대쪽) 행만 "
    #     "(같은 COHP# = 한 결합 블록, 예: Ru1–O19의 line 4–6, 8–26 …)"
    # )
    print(
        f"  {'COHP#':>5}  {'ele1':>4} {'idx1':>5}  "
        f"{'ele2':>4} {'idx2':>5}  {'distance':>8}  sum_ICOHP(sd×p)"
    )
    for cid in sorted(sums):
        if cid in meta:
            am0, an0, d0 = meta[cid]
            ele_mu, idx_mu = coarse_label_ele_idx(am0)
            ele_nu, idx_nu = coarse_label_ele_idx(an0)
            print(
                f"  {cid:5d}  {ele_mu:>4} {idx_mu:>5}  "
                f"{ele_nu:>4} {idx_nu:>5}  {d0:8.5f}  {sums[cid]:.5f}"
            )
        else:
            print(
                f"  {cid:5d}  {'?':>4} {'?':>5}  "
                f"{'?':>4} {'?':>5}  {'':>8}  {sums[cid]:.5f}"
            )

    if args.by_orbital_label:
        metal = sum_metal_sd_by_atom_mu(args.icohplist, args.spin)
        print()
        print("## 금속 atomMU (s, d) 라벨별 전역 합")
        print(f"# {'atomMU':<22}  sum_ICOHP")
        for k in sorted(metal, key=sort_key_mu):
            print(f"  {k:<22}  {metal[k]:.5f}")

    if args.by_oxygen_site:
        oxy = sum_oxygen_p_by_site(args.icohplist, args.spin)
        print()
        print("## 산소 site (p)")
        print(f"# {'site':<10}  sum_ICOHP")
        for k in sorted(oxy, key=sort_key_o):
            print(f"  {k:<10}  {oxy[k]:.5f}")


if __name__ == "__main__":
    main()
