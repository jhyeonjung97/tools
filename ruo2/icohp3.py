#!/usr/bin/env python3
"""LOBSTER ICOHPLIST.lobster 집계.

1) COHP#(같은 결합 블록)마다: 금속 s·d × 산소 p 행만 합산.
2) 금속 원자(Ru1, Ru2, …)마다: 산소 번호가 달라도 해당 금속이 끼인 모든 COHP# 합을 한 줄로 집계 (기본 출력).

옵션: 결합별 표, 라벨별·site별 집계.
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
    """Yields: cohp_id, atom_mu, atom_nu, distance, tx, ty, tz, icohp."""
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
                tx, ty, tz = int(parts[4]), int(parts[5]), int(parts[6])
                icohp = float(parts[7])
            except ValueError:
                continue
            yield cohp_id, parts[1], parts[2], dist, tx, ty, tz, icohp


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
) -> tuple[dict[int, float], dict[int, tuple[str, str, float, int, int, int]]]:
    """COHP# → (sd 금속 & p 산소) 행만 합산. 메타는 coarse Ru1–O19 행에서."""
    sums: dict[int, float] = defaultdict(float)
    meta: dict[int, tuple[str, str, float, int, int, int]] = {}

    for cohp_id, am, an, dist, tx, ty, tz, icohp in iter_icohp_rows(
        path, spin_filter
    ):
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
                meta[cohp_id] = (am, an, dist, tx, ty, tz)
            if (
                m_mu
                and m_nu
                and m_nu.group(1) != "O"
                and m_mu.group(1) == "O"
                and cohp_id not in meta
            ):
                meta[cohp_id] = (am, an, dist, tx, ty, tz)

        if row_metal_sd_and_oxygen_p(am, an):
            sums[cohp_id] += icohp

    return dict(sums), meta


def coarse_metal_label(meta_row: tuple[str, str, float, int, int, int]) -> str | None:
    """coarse (금속, 산소) 행에서 금속 쪽 라벨 Ru1, Hf16 등."""
    am0, an0 = meta_row[0], meta_row[1]
    m_mu, m_nu = ATOM_MU_RE.match(am0), ATOM_MU_RE.match(an0)
    if not m_mu or not m_nu:
        return None
    if m_mu.group(1) != "O" and m_nu.group(1) == "O":
        return am0
    if m_nu.group(1) != "O" and m_mu.group(1) == "O":
        return an0
    return None


def sum_by_metal_site(
    sums: dict[int, float], meta: dict[int, tuple[str, str, float, int, int, int]]
) -> dict[str, float]:
    """각 금속 원자에 대해 관련 COHP#(이웃 O가 달라도) sd×p 합."""
    out: dict[str, float] = defaultdict(float)
    for cid, val in sums.items():
        if cid not in meta:
            continue
        mlab = coarse_metal_label(meta[cid])
        if mlab:
            out[mlab] += val
    return dict(out)


def sum_metal_sd_by_atom_mu(path: Path, spin_filter: int | None) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    for _id, am, _an, *_rest, icohp in iter_icohp_rows(path, spin_filter):
        if is_metal_sd_atom_mu(am):
            sums[am] += icohp
    return dict(sums)


def sum_oxygen_p_by_site(path: Path, spin_filter: int | None) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    for _id, am, an, *_r, icohp in iter_icohp_rows(path, spin_filter):
        site = oxygen_site_from_row(am, an)
        if site:
            sums[site] += icohp
    return dict(sums)


def sort_key_mu(s: str) -> tuple:
    m = re.match(r"^([A-Za-z]+)(\d+)(?:_|$)", s)
    if m:
        el = m.group(1)
        pref = 0 if el == "Ru" else 1  # Ru 먼저, 그다음 Hf 등
        return (pref, el, int(m.group(2)), s)
    return (2, "", 0, s)


def sort_key_o(site: str) -> tuple:
    m = re.match(r"^O(\d+)$", site)
    if m:
        return (int(m.group(1)),)
    return (0,)


def main() -> None:
    p = argparse.ArgumentParser(
        description="ICOHPLIST: 금속 원자별 sd×p 합 (옵션: COHP#별, 라벨별, 산소 site)"
    )
    p.add_argument(
        "icohplist",
        type=Path,
        nargs="?",
        default=Path(__file__).resolve().parent / "ICOHPLIST.lobster",
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
        "--by-cohp",
        action="store_true",
        help="COHP#(결합)별 sd×p 합 표도 출력",
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
        for _id, am, _n, *_r, icohp in iter_icohp_rows(args.icohplist, args.spin):
            sums[am] += icohp
        print(f"# {args.icohplist}{spin_note}  —  atomMU별 sum(ICOH), 전 오비탈")
        print(f"# {'atomMU':<22}  sum_ICOHP")
        for k in sorted(sums, key=sort_key_mu):
            print(f"  {k:<22}  {sums[k]:.5f}")
        return

    sums, meta = sum_by_cohp_sd_mu_p_nu(args.icohplist, args.spin)
    by_metal = sum_by_metal_site(sums, meta)

    print(f"# {args.icohplist}{spin_note}")
    print(
        "# 금속 원자별: 각 COHP#에서 (금속 s·d × 산소 p)만 합한 뒤, "
        "같은 금속의 모든 이웃 결합을 더함 (산소 번호 무관)."
    )
    print(f"# {'metal':<8}  sum_ICOHP(sd×p, all neighbors)")
    for mlab in sorted(by_metal, key=sort_key_mu):
        print(f"  {mlab:<8}  {by_metal[mlab]:.5f}")

    if args.by_cohp:
        print()
        print(
            "## COHP#별 (한 결합씩): 금속 s·d × 산소 p "
            "(예: Ru1–O19 line 4–6, 8–26 …)"
        )
        print(
            f"# {'COHP#':>5}  {'atomMU':<6} {'atomNU':<6}  {'dist':>8}  "
            f"{'tx':>3}{'ty':>4}{'tz':>4}  sum_ICOHP(sd×p)"
        )
        for cid in sorted(sums):
            if cid in meta:
                am0, an0, d0, tx0, ty0, tz0 = meta[cid]
                print(
                    f"  {cid:5d}  {am0:<6} {an0:<6}  {d0:8.5f}  "
                    f"{tx0:3d}{ty0:4d}{tz0:4d}  {sums[cid]:.5f}"
                )
            else:
                print(
                    f"  {cid:5d}  {'?':<6} {'?':<6}  {'':>8}  {'':>11}  {sums[cid]:.5f}"
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
